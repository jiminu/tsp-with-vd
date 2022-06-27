#ifndef _FUNCTIONSFORVORONOIDIAGRAM3D_H
#define _FUNCTIONSFORVORONOIDIAGRAM3D_H

#include "rg_Const.h"
#include "Sphere.h"
#include "rg_Point3D.h"
#include "rg_Circle2D.h"
#include "rg_GeoFunc.h"

#include <math.h>


namespace V {

namespace GeometryTier {

const rg_INT XGE_YGE_ZGE = 1;
const rg_INT XGE_YGE_ZLT = 2;
const rg_INT XGE_YLT_ZGE = 3;
const rg_INT XGE_YLT_ZLT = 4;
const rg_INT XLT_YGE_ZGE = 5;
const rg_INT XLT_YGE_ZLT = 6;
const rg_INT XLT_YLT_ZGE = 7;
const rg_INT XLT_YLT_ZLT = 8;



template<class T> void sort3PointersInAscendingPowers(T** ptrs);
template<class T> void sort3PointersInAscendingPowers(T*& ptr1, T*& ptr2, T*& ptr3);
template<class T> void sort4PointersInAscendingPowers(T** ptrs);
template<class T> void sort4PointersInAscendingPowers(T*& ptr1, T*& ptr2, T*& ptr3, T*& ptr4);

rg_INT  classifyNormalSign(rg_REAL xNormal, rg_REAL yNormal, rg_REAL zNormal);
int compareInternalVoids(const void* ts1, const void* ts2);
int compareVDCellByID(const void* cell1, const void* cell2);
int compareBallsByRadius(const void* ball1, const void* ball2);
int compareBallsByRadiusInDescendingOrder(const void* ball1, const void* ball2);



rg_INT computeCircleTangentTo3CirclesOutside_GLOBAL( const rg_Circle2D& circle1, 
                                              const rg_Circle2D& circle2, 
                                              const rg_Circle2D& circle3, 
                                              rg_Circle2D* tangentCircle);
rg_INT computeMinTangentSpheresTo3SpheresOutside_GLOBAL( const Sphere& sphere1, 
                                                         const Sphere& sphere2, 
                                                         const Sphere& sphere3, 
                                                         Sphere* minTangentSphere);
Sphere computeMinSphereTangentToTwoBalls(const Sphere& ball1, const Sphere& ball2);



} // namespace GeometryTier

} // namespace V





template<class T> 
void V::GeometryTier::sort3PointersInAscendingPowers(T** ptrs)
{
//    unsigned int intTypePtrs[3] = { (unsigned int) ptrs[0],
//                                    (unsigned int) ptrs[1],
//                                    (unsigned int) ptrs[2] };
    T* intTypePtrs[3] = { ptrs[0], ptrs[1], ptrs[2] };

    if ( intTypePtrs[0] > intTypePtrs[1] )
        rg_SWAP(intTypePtrs[0], intTypePtrs[1]);

    if ( intTypePtrs[2] < intTypePtrs[0] )
    {
        rg_SWAP(intTypePtrs[1], intTypePtrs[2]);
        rg_SWAP(intTypePtrs[0], intTypePtrs[1]);
    }
    else if ( intTypePtrs[2] < intTypePtrs[1] )
    {
        rg_SWAP(intTypePtrs[1], intTypePtrs[2]);
    }
    
    for (int i=0; i<3; i++)
        ptrs[i] = (T*)intTypePtrs[i];
}

template<class T> 
void V::GeometryTier::sort3PointersInAscendingPowers(T*& ptr1, T*& ptr2, T*& ptr3)
{
    T* ptrArray[3] = { ptr1, ptr2, ptr3 };
    
    sort3PointersInAscendingPowers( ptrArray );

    ptr1 = ptrArray[0];
    ptr2 = ptrArray[1];
    ptr3 = ptrArray[2];
}


template<class T>
void V::GeometryTier::sort4PointersInAscendingPowers(T** ptrs)
{
    unsigned int intTypePtrs[4] = { (unsigned int) ptrs[0],
                                    (unsigned int) ptrs[1],
                                    (unsigned int) ptrs[2],
                                    (unsigned int) ptrs[3] };

    if ( intTypePtrs[0] > intTypePtrs[1] )
    {
        rg_SWAP(intTypePtrs[0], intTypePtrs[1]);
    }

    if ( intTypePtrs[2] > intTypePtrs[3] )
    {
        rg_SWAP(intTypePtrs[2], intTypePtrs[3]);
    }

    if ( intTypePtrs[0] < intTypePtrs[2] )
    {
        if ( intTypePtrs[1] > intTypePtrs[2] )
        {
            if ( intTypePtrs[1] < intTypePtrs[3] )
            {
                rg_SWAP(intTypePtrs[1], intTypePtrs[2]);                
            }
            else
            {
                rg_SWAP(intTypePtrs[1], intTypePtrs[2]);                
                rg_SWAP(intTypePtrs[2], intTypePtrs[3]);                
            }
        }
    }
    else 
    {
        if ( intTypePtrs[0] > intTypePtrs[3] )
        {
            rg_SWAP( intTypePtrs[0], intTypePtrs[2]);                
            rg_SWAP( intTypePtrs[1], intTypePtrs[3]);                
        }
        else
        {
            if ( intTypePtrs[1] < intTypePtrs[3] )
            {
                rg_SWAP( intTypePtrs[0], intTypePtrs[1]);                
                rg_SWAP( intTypePtrs[0], intTypePtrs[2]);                
            }
            else
            {
                rg_SWAP( intTypePtrs[0], intTypePtrs[1]);                
                rg_SWAP( intTypePtrs[0], intTypePtrs[2]);                
                rg_SWAP( intTypePtrs[2], intTypePtrs[3]);                
            }
        }
    }

    for (int i=0; i<4; i++)
    {
        ptrs[i] = (T*)intTypePtrs[i];
    }
}

template<class T>
void V::GeometryTier::sort4PointersInAscendingPowers(T*& ptr1, T*& ptr2, T*& ptr3, T*& ptr4)
{
    T* ptrArray[4] = { ptr1, ptr2, ptr3, ptr4 };
    
    sort4PointersInAscendingPowers( ptrArray );

    ptr1 = ptrArray[0];
    ptr2 = ptrArray[1];
    ptr3 = ptrArray[2];
    ptr4 = ptrArray[3];

}

#endif


