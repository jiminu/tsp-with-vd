#include "VIDICItem.h"
using namespace V::GeometryTier;


///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
VIDICItem::VIDICItem()
: m_vertex( rg_NULL ), m_AvailableNumVertexByBallConfiguration(0)
{
    for (rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        m_balls[i] = rg_NULL;
}




VIDICItem::VIDICItem( BallGenerator** balls, VDVertex* vertex )
: m_vertex( vertex ), m_AvailableNumVertexByBallConfiguration(0)
{
    for (rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        m_balls[i] = balls[i];
}


VIDICItem::VIDICItem( BallGenerator** balls, VDVertex* vertex, const rg_FLAG& numVertexBy4Balls )
: m_vertex( vertex ), m_AvailableNumVertexByBallConfiguration(numVertexBy4Balls)
{
    for (rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        m_balls[i] = balls[i];
}



VIDICItem::VIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4, VDVertex* vertex )
: m_vertex( vertex ), m_AvailableNumVertexByBallConfiguration(0)
{
    m_balls[0] = ball1;
    m_balls[1] = ball2;
    m_balls[2] = ball3;
    m_balls[3] = ball4;
}

VIDICItem::VIDICItem( const VIDICItem& anItem )
: m_vertex( anItem.m_vertex ), m_AvailableNumVertexByBallConfiguration(anItem.m_AvailableNumVertexByBallConfiguration)
{
    for (rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        m_balls[i] = anItem.m_balls[i];
}

VIDICItem::~VIDICItem()
{
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
BallGenerator* VIDICItem::getBall(const rg_INT& i) const
{
    if ( (i>-1) && (i<NUM_DEFINING_CELLS_OF_ONE_VERTEX) ) 
        return m_balls[i];
    else 
        return rg_NULL;

}

BallGenerator** VIDICItem::getAllBalls() 
{
    return m_balls;
}


VDVertex* VIDICItem::getVertex() const
{
    return m_vertex;
}

rg_FLAG VIDICItem::isSameVertexConfiguration(       BallGenerator* ball1, 
                                                    BallGenerator* ball2, 
                                                    BallGenerator* ball3, 
                                                    BallGenerator* ball4, 
                                              const rg_Point3D& vertexCoord )
{
    if (    ( m_balls[0] == ball1 ) && ( m_balls[1] == ball2 ) 
         && ( m_balls[2] == ball3 ) && ( m_balls[3] == ball4 ) )
    {
        if ( m_AvailableNumVertexByBallConfiguration == 1 )
            return rg_TRUE;
        else {        
            if ( m_vertex->getPoint() == vertexCoord )
                return rg_TRUE;
            else
                return rg_FALSE;
        }
    }
    else
    {
        return rg_FALSE;
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void VIDICItem::setBall( const rg_INT& i, BallGenerator* aBall )
{
    m_balls[i] = aBall;
}

void VIDICItem::setAllBalls( BallGenerator** balls )
{
    for (rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        m_balls[i] = balls[i];
}

rg_FLAG VIDICItem::setBall( BallGenerator* aBall )
{
    for (rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
    {
        if ( m_balls[i] == rg_NULL )
        {
            m_balls[i] = aBall;
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}


void VIDICItem::setVertex( VDVertex* aVertex )
{
    m_vertex = aVertex;
}


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
VIDICItem& VIDICItem::operator =( const VIDICItem& anItem )
{
    if ( this == &anItem )
        return *this;

    for (rg_INT i=0; i<NUM_DEFINING_CELLS_OF_ONE_VERTEX; i++)
        m_balls[i] = anItem.m_balls[i];

    m_vertex = anItem.m_vertex;
    m_AvailableNumVertexByBallConfiguration = anItem.m_AvailableNumVertexByBallConfiguration;

    return *this;
}

