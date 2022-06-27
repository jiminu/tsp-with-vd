#include "Polygon3D.h"



Polygon3D::Polygon3D()
{
}



Polygon3D::Polygon3D(const Polygon3D& polygon)
{
    duplicate(polygon.m_boundary_vertices);
}



Polygon3D::~Polygon3D()
{
}


void    Polygon3D::clear()
{
    m_boundary_vertices.clear();
}



Polygon3D& Polygon3D::operator =(const Polygon3D& polygon)
{
    if (this != &polygon) {
        duplicate(polygon.m_boundary_vertices);
    }

    return *this;
}



list< list<rg_Point3D> >& Polygon3D::boundary_vertices()
{
    return m_boundary_vertices;
}



const list< list<rg_Point3D> >& Polygon3D::boundary_vertices() const
{
    return m_boundary_vertices;
}



void Polygon3D::set_boundary_vertices(const list< list<rg_Point3D> >& boundaryVertices)
{
    duplicate(boundaryVertices);
}


void Polygon3D::duplicate(const list< list<rg_Point3D> >& boundaryVertices)
{
    m_boundary_vertices.clear();

    list< list<rg_Point3D> >::const_iterator i_list;
    for (i_list = boundaryVertices.begin(); i_list != boundaryVertices.end(); ++i_list) {
        m_boundary_vertices.push_back(*i_list);
    }
}

