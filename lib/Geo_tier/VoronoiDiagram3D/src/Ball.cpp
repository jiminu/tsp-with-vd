#include "Ball.h"
using namespace V::GeometryTier;

Ball::Ball()
: m_property(rg_NULL), m_IDFromInput(-1)
{
}



Ball::Ball(const rg_INT& ID)
: TopologicalEntity(ID), m_property(rg_NULL), m_IDFromInput(-1)
{
}



Ball::Ball(const Sphere& sphere, void* property, const rg_INT& IDFromInput)
: m_geometry(sphere), m_property(property), m_IDFromInput(IDFromInput)
{

}


Ball::Ball(const rg_INT& ID, const Sphere& sphere, void* property, const rg_INT& IDFromInput)
: TopologicalEntity(ID), m_geometry(sphere), m_property(property), m_IDFromInput(IDFromInput)
{
}



Ball::Ball(const Ball& ball)
: TopologicalEntity(ball), 
  m_geometry(ball.m_geometry), 
  m_property(ball.m_property), 
  m_IDFromInput(ball.m_IDFromInput)
{
}



Ball::~Ball()
{
    m_property = rg_NULL;
}




Sphere Ball::getGeometry() const
{
    return m_geometry;
}



void*  Ball::getProperty() const
{
    return m_property;
}



rg_INT Ball::getIDFromInput() const
{
    return m_IDFromInput;
}



void   Ball::setGeometry(const Sphere& sphere)
{
    m_geometry = sphere;
}



void   Ball::setProperty(void* property)
{
    m_property = property;
}



void Ball::setIDFromInput(const rg_INT& IDFromInput)
{
    m_IDFromInput = IDFromInput;
}



Ball&  Ball::operator=(const Ball& ball)
{
    if ( this == &ball ) {
        return *this;
    }

    TopologicalEntity::operator =(ball);

    m_geometry = ball.m_geometry;
    m_property = ball.m_property;
    m_IDFromInput = ball.m_IDFromInput;

    return *this;
}



