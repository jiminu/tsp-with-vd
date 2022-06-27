#ifndef _VIDIC_H
#define _VIDIC_H

#include "rg_dList.h"

#include "rg_BallGenerator.h"
#include "VDVertex.h"

#include "VIDICItem.h"

namespace V {

namespace GeometryTier {



const rg_REAL THRESHOLD_FOR_REHASHING = 0.9;

class VIDIC
{
private:
    rg_INT                 m_numOfItems;
    rg_INT                 m_numOfItemsForRehashing;

    rg_INT                 m_size;
    rg_dList< VIDICItem >* m_hashTable;

public:
    //  constructor & deconstructor..
    VIDIC();
    VIDIC( const rg_INT& size );
    VIDIC( const VIDIC& aVIDIC );
    ~VIDIC();

    //  get functions.. 
    rg_INT  getNumOfItems() const;
    rg_INT  getSize() const;

    rg_dList< VIDICItem >* getHashTable();

    //  set functions..
    void setSize( const rg_INT& size );

    //  operator overloading..
    VIDIC& operator =( const VIDIC& aVIDIC );


    //  functions for hash table
    VIDICItem* findVIDICItem( const VIDICItem& anItem );
    VIDICItem* findVIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, VDVertex* vertex );
    VIDICItem* findVIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, const rg_Point3D& vertexCoord );

    VIDICItem* insertVIDICItem( const VIDICItem& anItem );
    VIDICItem* insertVIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, VDVertex* vertex );
    VIDICItem* insertVIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, VDVertex* vertex, const rg_FLAG& numVertexBy4Balls );

    rg_FLAG    removeVIDICItem( const VIDICItem& anItem );
    rg_FLAG    removeVIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, VDVertex* vertex );

private:
    //  private functions for hash table
    void   sort4BallPointers( BallGenerator** ptrBalls );
    void   sort4BallPointers( BallGenerator*& ptrBall1, BallGenerator*& ptrBall2, BallGenerator*& ptrBall3, BallGenerator*& ptrBall4 );
        //void swap( BallGenerator*&      ptrBall1,    BallGenerator*&      ptrBall2,
        //           unsigned int& intPtrBall1, unsigned int& intPtrBall2 );
        void swap( BallGenerator*&      ptrBall1,    BallGenerator*&      ptrBall2,
                   size_t& intPtrBall1, size_t& intPtrBall2 );
    rg_INT evaluateHashFunction( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4 );
    void   rehashVIDICByDoubleSize();
    void   removeAll();

};


} // namespace GeometryTier

} // namespace V


#endif
