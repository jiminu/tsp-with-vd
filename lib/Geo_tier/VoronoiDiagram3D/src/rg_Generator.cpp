#include "rg_Generator.h"
#include "VDCell.h"
using namespace V::GeometryTier;


///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
Generator::Generator()
: m_ID(-1), m_cell(rg_NULL), m_checkForBucket(rg_FALSE)
{
}

Generator::Generator(const rg_INT& ID)
: m_ID(ID), m_cell(rg_NULL), m_checkForBucket(rg_FALSE)
{
}

Generator::Generator(const rg_INT& ID, VDCell* cell)
: m_ID(ID), m_cell(cell), m_checkForBucket(rg_FALSE)
{
}

Generator::Generator(const Generator& aGenerator)
: m_ID(aGenerator.m_ID), m_cell(aGenerator.m_cell), m_checkForBucket(aGenerator.m_checkForBucket)
{
}

Generator::~Generator()
{
}


///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
rg_INT  Generator::getID() const
{
    return m_ID;
}

VDCell* Generator::getCell() const
{
    return m_cell;
}

rg_INT Generator::getTypeOfGenerator() const
{
    return NO_TYPE;
}
    
rg_FLAG Generator::isThisGeneratorIntersectedWith(Generator* anotherGenerator)
{
    return rg_UNKNOWN;
}

rg_FLAG Generator::isThisGeneratorIntersectedWith(const Sphere& aSphere)
{
    return rg_UNKNOWN;
}

Generator* Generator::duplicate() const
{
    Generator* aCopy = new Generator( *this );

    return aCopy;
}

rg_FLAG Generator::isTangible() const
{
    return m_ProbeTangibility;
}

void    Generator::isTangible(const rg_FLAG& tangibility)
{
    m_ProbeTangibility = tangibility;
}

///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void Generator::setID(const rg_INT& ID)
{
    m_ID = ID;
}

void Generator::setCell(VDCell* cell)
{
    m_cell = cell;
}

void Generator::setGenerator(const rg_INT& ID, VDCell* cell)
{
    m_ID   = ID;
    m_cell = cell;
}

void Generator::setGenerator(const Generator& aGenerator)
{
    m_ID             = aGenerator.m_ID;
    m_cell           = aGenerator.m_cell;
    m_checkForBucket = aGenerator.m_checkForBucket;
}


///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
Generator& Generator::operator =(const Generator& aGenerator)
{
    if ( this == &aGenerator )
        return *this;

    m_ID             = aGenerator.m_ID;
    m_cell           = aGenerator.m_cell;
    m_checkForBucket = aGenerator.m_checkForBucket;

    return *this;
}

