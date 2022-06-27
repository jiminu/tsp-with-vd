#include "VertexQT2D.h"
using namespace V::GeometryTier;


VertexQT2D::VertexQT2D()
: m_disc(rg_NULL)
{
}



VertexQT2D::VertexQT2D(const rg_INT& ID, Disc* disc)
: VertexSCDS(ID),
  m_disc(disc)
{
}



VertexQT2D::VertexQT2D(const VertexQT2D& vertex)
: VertexSCDS(vertex),
  m_disc(vertex.m_disc)
{
}



VertexQT2D::~VertexQT2D()
{
}




Disc* VertexQT2D::getDisc()
{
    return m_disc;
}



rg_Circle2D VertexQT2D::getDiscGeometry() const
{
    return m_disc->getGeometry();
}



rg_BOOL VertexQT2D::isVirtual() const
{
    if ( m_disc == rg_NULL ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



void VertexQT2D::setDisc(Disc* disc)
{
    m_disc = disc;
}




VertexQT2D& VertexQT2D::operator =(const VertexQT2D& vertex)
{
    if ( this == &vertex ) {
        return *this;
    }

    VertexSCDS::operator =(vertex);

    m_disc = vertex.m_disc;

    return *this;
}



