#include "rg_BallGenerator.h"
using namespace V::GeometryTier;

///////////////////////////////////////////////////////////////////////////////
//
//  constructor & deconstructor..
BallGenerator::BallGenerator()
: m_ID(-1), m_property( rg_NULL ), m_IDFromInput(-1), m_checkForBucket(rg_FALSE)
{
}



BallGenerator::BallGenerator(const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, const rg_REAL& radius, void* property, const rg_INT& IDFromInput)
: m_ID(-1), Sphere( x, y, z, radius ), m_property( property ), m_IDFromInput(IDFromInput), m_checkForBucket(rg_FALSE)
{
}



BallGenerator::BallGenerator(const rg_Point3D& center, const rg_REAL& radius, void*  property, const rg_INT& IDFromInput)
: m_ID(-1), Sphere( center, radius ), m_property( property ), m_IDFromInput(IDFromInput), m_checkForBucket(rg_FALSE)
{
}



BallGenerator::BallGenerator(const rg_INT& ID, const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, const rg_REAL& radius, void* property, const rg_INT& IDFromInput)
: m_ID(ID), Sphere( x, y, z, radius ), m_property( property ), m_IDFromInput(IDFromInput), m_checkForBucket(rg_FALSE)
{
}



BallGenerator::BallGenerator(const rg_INT& ID, const rg_Point3D& center, const rg_REAL& radius, void* property, const rg_INT& IDFromInput)
: m_ID(ID), Sphere( center, radius ), m_property( property ), m_IDFromInput(IDFromInput), m_checkForBucket(rg_FALSE)
{
}



BallGenerator::BallGenerator(const rg_INT& ID,  const Sphere& ball, void* property, const rg_INT& IDFromInput)
: m_ID(ID), Sphere( ball ), m_property( property ), m_IDFromInput(IDFromInput), m_checkForBucket(rg_FALSE)
{
}



BallGenerator::BallGenerator(const BallGenerator& aBallGenerator)
{
    m_ID        = aBallGenerator.m_ID;
    m_cell      = aBallGenerator.m_cell;

    m_center    = aBallGenerator.m_center;
    m_radius    = aBallGenerator.m_radius;
    
    m_property  = aBallGenerator.m_property;
    m_IDFromInput = aBallGenerator.m_IDFromInput;

    m_checkForBucket = aBallGenerator.m_checkForBucket;
}



BallGenerator::~BallGenerator()
{
}



///////////////////////////////////////////////////////////////////////////////
//
//  get functions.. 
rg_INT  BallGenerator::getID() const
{
    return m_ID;
}



VDCell* BallGenerator::getCell() const
{
    return m_cell;
}



Sphere BallGenerator::getBall() const
{
    return Sphere(m_center, m_radius);
}



void* BallGenerator::getProperty() const
{
    return m_property;
}



rg_INT BallGenerator::getIDFromInput() const
{
    return m_IDFromInput;
}



rg_FLAG BallGenerator::isThisGeneratorIntersectedWith(BallGenerator* anotherGenerator)
{
    return this->isThereIntersectionWith( *anotherGenerator );
}



rg_FLAG BallGenerator::isThisGeneratorIntersectedWith(const Sphere& aSphere)
{
    return this->isThereIntersectionWith( aSphere );
}


/*
Generator* BallGenerator::duplicate() const
{
    BallGenerator* aCopy = new BallGenerator( *this );

    return (Generator*) aCopy;
}
*/


rg_FLAG BallGenerator::isTangible() const
{
    return m_ProbeTangibility;
}

void    BallGenerator::isTangible(const rg_FLAG& tangibility)
{
    m_ProbeTangibility = tangibility;
}


///////////////////////////////////////////////////////////////////////////////
//
//  set functions..
void BallGenerator::setID(const rg_INT& ID)
{
    m_ID = ID;
}



void BallGenerator::setCell(VDCell* cell)
{
    m_cell = cell;
}



void BallGenerator::setProperty(void* property)
{
    m_property = property;
}



void BallGenerator::setIDFromInput(const rg_INT& IDFromInput)
{
    m_IDFromInput = IDFromInput;
}



///////////////////////////////////////////////////////////////////////////////
//
//  operator overloading..
BallGenerator& BallGenerator::operator =(const BallGenerator& aBallGenerator)
{
    if ( this == &aBallGenerator )
        return *this;

    m_ID        = aBallGenerator.m_ID;
    m_cell      = aBallGenerator.m_cell;

    m_center    = aBallGenerator.m_center;
    m_radius    = aBallGenerator.m_radius;
    
    m_property  = aBallGenerator.m_property;
    m_IDFromInput = aBallGenerator.m_IDFromInput;

    m_checkForBucket = aBallGenerator.m_checkForBucket;

    return *this;
}

