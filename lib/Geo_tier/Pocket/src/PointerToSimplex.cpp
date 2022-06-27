#include "PointerToSimplex.h"

PointerToSimplex::PointerToSimplex()
{
    m_typeOfSimplex  = rg_UNKNOWN;
    m_simplex        = rg_NULL;
}



PointerToSimplex::PointerToSimplex(const rg_INT& typeOfSimplex, void* simplex)
{
    m_typeOfSimplex = typeOfSimplex;
    m_simplex       = simplex;
}


PointerToSimplex::PointerToSimplex(const PointerToSimplex& tempPointerToSimplex)
{
    m_typeOfSimplex = tempPointerToSimplex.m_typeOfSimplex;
    m_simplex       = tempPointerToSimplex.m_simplex;
}


PointerToSimplex::PointerToSimplex(AugmentedBCVertex* vertex)
{
	m_typeOfSimplex = VERTEX_SIMPLEX;
	m_simplex = (void*)vertex;
}

PointerToSimplex::PointerToSimplex(AugmentedBCEdge* edge)
{
	m_typeOfSimplex = EDGE_SIMPLEX;
	m_simplex = (void*)edge;
}


PointerToSimplex::PointerToSimplex(AugmentedBCFace* face)
{
	m_typeOfSimplex = FACE_SIMPLEX;
	m_simplex = (void*)face;
}


PointerToSimplex::PointerToSimplex(AugmentedBCCell* cell)
{
	m_typeOfSimplex = CELL_SIMPLEX;
	m_simplex = (void*)cell;
}


PointerToSimplex::PointerToSimplex(BetaVertex* vertex)
{
	m_typeOfSimplex = VERTEX_SIMPLEX;
	m_simplex = (void*)vertex;
}


PointerToSimplex::PointerToSimplex(BetaEdge* edge)
{
	m_typeOfSimplex = EDGE_SIMPLEX;
	m_simplex = (void*)edge;
}


PointerToSimplex::PointerToSimplex(BetaFace* face)
{
	m_typeOfSimplex = FACE_SIMPLEX;
	m_simplex = (void*)face;
}



PointerToSimplex::PointerToSimplex(BetaCell* cell)
{
	m_typeOfSimplex = CELL_SIMPLEX;
	m_simplex = (void*)cell;
}


PointerToSimplex::~PointerToSimplex()
{
}

 
rg_INT PointerToSimplex::getTypeOfSimplex() const
{
    return m_typeOfSimplex;
}


void* PointerToSimplex::getSimplex() const
{
    return m_simplex;
}



void  PointerToSimplex::setTypeOfSimplex(const rg_INT& typeOfSimplex)
{
    m_typeOfSimplex = typeOfSimplex;
}



void PointerToSimplex::setSimplex(void* simplex)
{
    m_simplex       = simplex;
}




void PointerToSimplex::setPointerToSimplex(const rg_INT& typeOfSimplex, void* simplex)
{
     m_typeOfSimplex = typeOfSimplex;
     m_simplex       = simplex;
}

PointerToSimplex& PointerToSimplex::operator =(const PointerToSimplex& tempPointerToSimplex)
{
    if ( this == &tempPointerToSimplex )
        return *this;

    m_typeOfSimplex = tempPointerToSimplex.m_typeOfSimplex;
    m_simplex       = tempPointerToSimplex.m_simplex;

    return *this;

}


PointerToSimplex& PointerToSimplex::operator =(AugmentedBCVertex* vertex)
{
    m_typeOfSimplex = VERTEX_SIMPLEX;
    m_simplex       = (void*)vertex;

    return *this;
}



PointerToSimplex& PointerToSimplex::operator =(AugmentedBCEdge* edge)
{
    m_typeOfSimplex = EDGE_SIMPLEX;
    m_simplex       = (void*)edge;

    return *this;
}



PointerToSimplex& PointerToSimplex::operator =(AugmentedBCFace* face)
{
    m_typeOfSimplex = FACE_SIMPLEX;
    m_simplex       = (void*)face;

    return *this;
}


PointerToSimplex& PointerToSimplex::operator =(AugmentedBCCell* cell)
{
    m_typeOfSimplex = CELL_SIMPLEX;
    m_simplex       = (void*)cell;

    return *this;
}
