#include "TopologicalEntity.h"

///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
TopologicalEntity::TopologicalEntity()
: m_ID(-1)
, m_delete(false)
, m_flag_query( false )
, m_flag_query_TOI( false )
{
}

TopologicalEntity::TopologicalEntity(const rg_INT& ID)
: m_ID(ID)
, m_delete(false)
, m_flag_query( false )
, m_flag_query_TOI( false )
{
}

TopologicalEntity::TopologicalEntity(const TopologicalEntity& aTopoEntity)
: m_ID(aTopoEntity.m_ID)
, m_delete(aTopoEntity.m_delete)
, m_flag_query( aTopoEntity.m_flag_query )
, m_flag_query_TOI( aTopoEntity.m_flag_query_TOI )
{
}

TopologicalEntity::~TopologicalEntity()
{
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
//rg_INT  TopologicalEntity::getID() const
//{
//    return m_ID;
//}


///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void TopologicalEntity::setID(const rg_INT& ID)
{
    m_ID = ID;
}

void TopologicalEntity::setTopologicalEntity(const TopologicalEntity& aTopoEntity)
{
    m_ID = aTopoEntity.m_ID;
}


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
TopologicalEntity& TopologicalEntity::operator =(const TopologicalEntity& aTopoEntity)
{
    if ( this == &aTopoEntity )
        return *this;

    m_ID = aTopoEntity.m_ID;
    m_delete = aTopoEntity.m_delete;
    m_flag_query = aTopoEntity.m_flag_query;
    m_flag_query_TOI = aTopoEntity.m_flag_query_TOI;

    return *this;
}

