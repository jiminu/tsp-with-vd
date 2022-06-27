#ifndef POLYGON3D_H
#define POLYGON3D_H

#include "rg_Point3D.h"

#include <list>
using namespace std;


class Polygon3D
{
private:
    list< list<rg_Point3D> > m_boundary_vertices;


public:
    Polygon3D();
    Polygon3D(const Polygon3D& polygon);
    ~Polygon3D();

    void    clear();

    Polygon3D& operator =(const Polygon3D& polygon);

    list< list<rg_Point3D> >& boundary_vertices();

    const list< list<rg_Point3D> >& boundary_vertices() const;

    void set_boundary_vertices(const list< list<rg_Point3D> >& boundaryVertices);

private:
    void duplicate(const list< list<rg_Point3D> >& boundaryVertices);
};

#endif


