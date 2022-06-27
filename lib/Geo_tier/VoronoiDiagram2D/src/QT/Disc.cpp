#include "Disc.h"
using namespace V::GeometryTier;


Disc::Disc()
: m_property(rg_NULL)
{
}



Disc::Disc(const rg_INT& ID)
: TopologicalEntity(ID),
  m_property(rg_NULL)
{
}



Disc::Disc(const rg_INT& ID, const rg_Circle2D& circle, void* discProperty)
: TopologicalEntity(ID),
  m_circle(circle),
  m_property(discProperty)
{
}



Disc::Disc(const Disc& disc)
: TopologicalEntity(disc),
  m_circle(disc.m_circle),
  m_property(disc.m_property)
{
}



Disc::~Disc()
{
}



rg_Circle2D Disc::getGeometry() const
{
    return m_circle;
}



void* Disc::getProperty() const
{
    return m_property;
}




void Disc::setGeometry(const rg_Circle2D& circle)
{
    m_circle = circle;
}



void Disc::setProperty(void* discProperty)
{
    m_property = discProperty;
}




Disc& Disc::operator=(const Disc& disc)
{
    if ( this == &disc ) {
        return *this;
    }

    TopologicalEntity::operator =(disc);
    m_circle   = disc.m_circle;
    m_property = disc.m_property;

    return *this;
}



