#include "VIDIC.h"
using namespace V::GeometryTier;



///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VIDIC::VIDIC()
: m_numOfItems(0), m_numOfItemsForRehashing(0), m_size(0), m_hashTable(rg_NULL)
{
}

VIDIC::VIDIC( const rg_INT& size )
: m_numOfItems(0), m_size(size)
{
    m_numOfItemsForRehashing = (rg_INT)(THRESHOLD_FOR_REHASHING * m_size);

    m_hashTable              = new rg_dList< VIDICItem >[ m_size ];
}

VIDIC::VIDIC( const VIDIC& aVIDIC )
: m_numOfItems( aVIDIC.m_numOfItems ), 
  m_numOfItemsForRehashing( aVIDIC.m_numOfItemsForRehashing ), 
  m_size( aVIDIC.m_size )
{
    m_hashTable = new rg_dList< VIDICItem >[ m_size ];

    for (rg_INT i=0; i<m_size; i++)
        m_hashTable[i].duplicateList( aVIDIC.m_hashTable[i] );

}

VIDIC::~VIDIC()
{
    removeAll();
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
rg_INT VIDIC::getNumOfItems() const
{
    return m_numOfItems;
}

rg_INT VIDIC::getSize() const
{
    return m_size;
}


rg_dList< VIDICItem >* VIDIC::getHashTable()
{
    return m_hashTable;
}


///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VIDIC::setSize( const rg_INT& size )
{
    removeAll();

    m_size                   = size;
    m_numOfItemsForRehashing = (rg_INT) (THRESHOLD_FOR_REHASHING * m_size);
    m_hashTable              = new rg_dList< VIDICItem >[ m_size ];
}


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
VIDIC& VIDIC::operator =( const VIDIC& aVIDIC )
{
    if ( this == &aVIDIC )
        return *this;

    removeAll();

    m_numOfItems             = aVIDIC.m_numOfItems;
    m_numOfItemsForRehashing = aVIDIC.m_numOfItemsForRehashing;
    m_size                   = aVIDIC.m_size;
    m_hashTable              = new rg_dList< VIDICItem >[ m_size ];

    for (rg_INT i=0; i<m_size; i++)
        m_hashTable[i].duplicateList( aVIDIC.m_hashTable[i] );

    return *this;
}



///////////////////////////////////////////////////////////////////////////////
//
//  functions for hash table
VIDICItem* VIDIC::findVIDICItem( const VIDICItem& anItem )
{
    VIDICItem* theItem = findVIDICItem( anItem.m_balls[0], 
                                        anItem.m_balls[1], 
                                        anItem.m_balls[2], 
                                        anItem.m_balls[3], 
                                        anItem.m_vertex->getPoint() );

    return theItem;
}


VIDICItem* VIDIC::findVIDICItem( BallGenerator*   ball1, 
                                 BallGenerator*   ball2, 
                                 BallGenerator*   ball3, 
                                 BallGenerator*   ball4, 
                                 VDVertex* vertex)
{
    VIDICItem* theItem = findVIDICItem( ball1, ball2, ball3, ball4, vertex->getPoint() );

    return theItem;
}

VIDICItem* VIDIC::findVIDICItem(       BallGenerator*     ball1, 
                                       BallGenerator*     ball2, 
                                       BallGenerator*     ball3, 
                                       BallGenerator*     ball4, 
                                 const rg_Point3D& vertexCoord )
{
    BallGenerator* ptrBalls[NUM_DEFINING_CELLS_OF_ONE_VERTEX] = { ball1, ball2, ball3, ball4 };
    
    sort4BallPointers( ptrBalls );

    rg_INT tableIndex = evaluateHashFunction( ptrBalls[0], ptrBalls[1], ptrBalls[2], ptrBalls[3] );

    VIDICItem* currentItem = rg_NULL;

	m_hashTable[ tableIndex ].reset4Loop();
	while( m_hashTable[ tableIndex ].setNext4Loop() )
	{
		currentItem = m_hashTable[ tableIndex ].getpEntity();

        if ( currentItem->isSameVertexConfiguration( ptrBalls[0], ptrBalls[1], ptrBalls[2], ptrBalls[3], vertexCoord) )
        {
            return currentItem;
        }
    }
    
    return rg_NULL;
}

VIDICItem* VIDIC::insertVIDICItem( const VIDICItem& anItem )
{
    VIDICItem* ptrInsertingItem =  insertVIDICItem( anItem.m_balls[0], 
                                                    anItem.m_balls[1], 
                                                    anItem.m_balls[2], 
                                                    anItem.m_balls[3], 
                                                    anItem.m_vertex   ); 
    return ptrInsertingItem;
}

VIDICItem* VIDIC::insertVIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, VDVertex* vertex)
{
    BallGenerator* ptrBalls[NUM_DEFINING_CELLS_OF_ONE_VERTEX] = { ball1, ball2, ball3, ball4 };
    
    sort4BallPointers( ptrBalls );

    rg_INT tableIndex = evaluateHashFunction( ptrBalls[0], ptrBalls[1], ptrBalls[2], ptrBalls[3] );

    VIDICItem itemToInsert( ptrBalls, vertex );

    VIDICItem* ptrInsertingItem = m_hashTable[ tableIndex ].addTail( itemToInsert );

    m_numOfItems++;

    if ( m_numOfItems > m_numOfItemsForRehashing )
    {
        rehashVIDICByDoubleSize();
    }

    return ptrInsertingItem;
}



VIDICItem* VIDIC::insertVIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, VDVertex* vertex, const rg_FLAG& numVertexBy4Balls )
{
    BallGenerator* ptrBalls[NUM_DEFINING_CELLS_OF_ONE_VERTEX] = { ball1, ball2, ball3, ball4 };
    
    sort4BallPointers( ptrBalls );

    rg_INT tableIndex = evaluateHashFunction( ptrBalls[0], ptrBalls[1], ptrBalls[2], ptrBalls[3] );

    VIDICItem itemToInsert( ptrBalls, vertex, numVertexBy4Balls );

    VIDICItem* ptrInsertingItem = m_hashTable[ tableIndex ].addTail( itemToInsert );

    m_numOfItems++;

    if ( m_numOfItems > m_numOfItemsForRehashing )
    {
        rehashVIDICByDoubleSize();
    }

    return ptrInsertingItem;
}




rg_FLAG VIDIC::removeVIDICItem( const VIDICItem& anItem )
{
    rg_FLAG isRemoved =  removeVIDICItem( anItem.m_balls[0], 
                                          anItem.m_balls[1], 
                                          anItem.m_balls[2], 
                                          anItem.m_balls[3], 
                                          anItem.m_vertex   ); 
    return isRemoved;
}

rg_FLAG VIDIC::removeVIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, VDVertex* vertex )
{
    BallGenerator* ptrBalls[NUM_DEFINING_CELLS_OF_ONE_VERTEX] = { ball1, ball2, ball3, ball4 };
    
    sort4BallPointers( ptrBalls );

    rg_INT tableIndex = evaluateHashFunction( ptrBalls[0], ptrBalls[1], ptrBalls[2], ptrBalls[3] );



    rg_INT sizeOfHashTable = m_hashTable[ tableIndex ].getSize();

    rg_dNode< VIDICItem >* currNode = m_hashTable[ tableIndex ].getHead();
    rg_dNode< VIDICItem >* nextNode = rg_NULL;

    for ( rg_INT i=0; i<sizeOfHashTable; i++)
    {
        nextNode = currNode->getNext();
       
        if ( currNode->getpEntity()->isSameVertexConfiguration( ptrBalls[0], ptrBalls[1], ptrBalls[2], ptrBalls[3], vertex->getPoint() ) )
        {
            m_hashTable[ tableIndex ].kill( currNode );
            m_numOfItems--;

            return rg_TRUE;
        }
        currNode = nextNode;
    }

/*
    VIDICItem* currentItem = rg_NULL;

	m_hashTable[ tableIndex ].reset4Loop();
	while( m_hashTable[ tableIndex ].setNext4Loop() )
	{
		currentItem = m_hashTable[ tableIndex ].getpEntity();

        if ( currentItem->isSameVertexConfiguration( ptrBalls[0], ptrBalls[1], ptrBalls[2], ptrBalls[3], vertex->getPoint() ) )
        {
            m_hashTable[ tableIndex ].killCurrent();

            m_numOfItems--;

            return rg_TRUE;
        }
    }
*/
    return rg_FALSE;
}

///////////////////////////////////////////////////////////////////////////////
//
//  private functions for hash table
void VIDIC::sort4BallPointers( BallGenerator** ptrBalls )
{
    size_t intPtrBall[NUM_DEFINING_CELLS_OF_ONE_VERTEX] 
                     = { (size_t) ptrBalls[0],
                         (size_t) ptrBalls[1],
                         (size_t) ptrBalls[2],
                         (size_t) ptrBalls[3] };

    if ( intPtrBall[0] > intPtrBall[1] )
    {
        swap(ptrBalls[0], ptrBalls[1], intPtrBall[0], intPtrBall[1]);
    }

    if ( intPtrBall[2] > intPtrBall[3] )
    {
        swap(ptrBalls[2], ptrBalls[3], intPtrBall[2], intPtrBall[3]);
    }

    if ( intPtrBall[0] < intPtrBall[2] )
    {
        if ( intPtrBall[1] > intPtrBall[2] )
        {
            if ( intPtrBall[1] < intPtrBall[3] )
            {
                swap(ptrBalls[1], ptrBalls[2], intPtrBall[1], intPtrBall[2]);                
            }
            else
            {
                swap(ptrBalls[1], ptrBalls[2], intPtrBall[1], intPtrBall[2]);                
                swap(ptrBalls[2], ptrBalls[3], intPtrBall[2], intPtrBall[3]);                
            }
        }
    }
    else 
    {
        if ( intPtrBall[0] > intPtrBall[3] )
        {
            swap(ptrBalls[0], ptrBalls[2], intPtrBall[0], intPtrBall[2]);                
            swap(ptrBalls[1], ptrBalls[3], intPtrBall[1], intPtrBall[3]);                
        }
        else
        {
            if ( intPtrBall[1] < intPtrBall[3] )
            {
                swap(ptrBalls[0], ptrBalls[1], intPtrBall[0], intPtrBall[1]);                
                swap(ptrBalls[0], ptrBalls[2], intPtrBall[0], intPtrBall[2]);                
            }
            else
            {
                swap(ptrBalls[0], ptrBalls[1], intPtrBall[0], intPtrBall[1]);                
                swap(ptrBalls[0], ptrBalls[2], intPtrBall[0], intPtrBall[2]);                
                swap(ptrBalls[2], ptrBalls[3], intPtrBall[2], intPtrBall[3]);                
            }
        }
    }
}

void VIDIC::sort4BallPointers( BallGenerator*& ptrBall1, BallGenerator*& ptrBall2, BallGenerator*& ptrBall3, BallGenerator*& ptrBall4 )
{
    BallGenerator* ptrBalls[NUM_DEFINING_CELLS_OF_ONE_VERTEX] = { ptrBall1, ptrBall2, ptrBall3, ptrBall4 };

    sort4BallPointers( ptrBalls );

    ptrBall1 = ptrBalls[0];
    ptrBall2 = ptrBalls[1];
    ptrBall3 = ptrBalls[2];
    ptrBall4 = ptrBalls[3];
}

//void VIDIC::swap( BallGenerator*&      ptrBall1,    BallGenerator*&      ptrBall2,
//                  unsigned int& intPtrBall1, unsigned int& intPtrBall2 )
//{
//    BallGenerator*      tempPtrBall    = ptrBall1;
//    unsigned int tempIntPtrBall = intPtrBall1;
//
//    ptrBall1 = ptrBall2;
//    ptrBall2 = tempPtrBall;
//
//    intPtrBall1 = intPtrBall2;
//    intPtrBall2 = tempIntPtrBall;
//}

void VIDIC::swap( BallGenerator*&      ptrBall1,    BallGenerator*&      ptrBall2,
                  size_t& intPtrBall1, size_t& intPtrBall2 )
{
    BallGenerator*      tempPtrBall    = ptrBall1;
    unsigned int tempIntPtrBall = intPtrBall1;

    ptrBall1 = ptrBall2;
    ptrBall2 = tempPtrBall;

    intPtrBall1 = intPtrBall2;
    intPtrBall2 = tempIntPtrBall;
}



rg_INT VIDIC::evaluateHashFunction( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4 )
{
    size_t hashcode = 0;
    
    size_t integerOfPtr[NUM_DEFINING_CELLS_OF_ONE_VERTEX]
                     = { (size_t) ball1,
                         (size_t) ball2,
                         (size_t) ball3,
                         (size_t) ball4 };
    for ( rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        hashcode += integerOfPtr[i];

    //  using "division method"
    rg_INT valueOfCompressionMap = hashcode%m_size;

    return valueOfCompressionMap;


    //unsigned long hashcode = 0;
    //
    //unsigned int integerOfPtr[NUM_DEFINING_CELLS_OF_ONE_VERTEX] 
    //                 = { (unsigned int) ball1, 
    //                     (unsigned int) ball2, 
    //                     (unsigned int) ball3, 
    //                     (unsigned int) ball4 };
    //for ( rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
    //    hashcode += integerOfPtr[i];

    ////  using "division method"
    //rg_INT valueOfCompressionMap = hashcode%m_size;

    //return valueOfCompressionMap;
}

void VIDIC::rehashVIDICByDoubleSize()
{
    rg_INT sizeBeforeRehashing = m_size;
    
    m_size                   = sizeBeforeRehashing*2;
    m_numOfItemsForRehashing = (rg_INT) (THRESHOLD_FOR_REHASHING * m_size);

    rg_dList< VIDICItem >* hashTableToBeRehashed = new rg_dList< VIDICItem >[ m_size ];

    for (rg_INT i=0; i<sizeBeforeRehashing; i++)
    {
        VIDICItem* currentItem = rg_NULL;

	    m_hashTable[ i ].reset4Loop();
	    while( m_hashTable[ i ].setNext4Loop() )
	    {
		    currentItem = m_hashTable[ i ].getpEntity();

            rg_INT tableIndex = evaluateHashFunction( currentItem->m_balls[0], 
                                                      currentItem->m_balls[1], 
                                                      currentItem->m_balls[2], 
                                                      currentItem->m_balls[3] ); 

            hashTableToBeRehashed[ tableIndex ].addTail( *currentItem );
        }
    }

    if ( m_hashTable != rg_NULL )
        delete [] m_hashTable;

    m_hashTable = hashTableToBeRehashed;
}



void VIDIC::removeAll()
{
    m_numOfItems = 0;
    m_size       = 0;

    if ( m_hashTable != rg_NULL )
        delete [] m_hashTable;
}

