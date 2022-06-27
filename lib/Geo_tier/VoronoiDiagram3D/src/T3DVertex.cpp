#include "T3DVertex.h"
using namespace V::GeometryTier;


T3DVertex::T3DVertex()
: m_firstTetrahedron(rg_NULL), m_property(rg_NULL), m_check(NOT_VISITED)
{
}



T3DVertex::T3DVertex(const rg_Point3D& point)
: m_point(point), m_firstTetrahedron(rg_NULL), m_property(rg_NULL), m_check(NOT_VISITED)
{
}



T3DVertex::T3DVertex(const rg_INT& ID, const rg_Point3D& point)
: TopologicalEntity(ID), m_point(point), m_firstTetrahedron(rg_NULL), m_property(rg_NULL), m_check(NOT_VISITED)
{
}



T3DVertex::T3DVertex(const T3DVertex& vertex)
{
	m_ID               = vertex.m_ID;
	m_point            = vertex.m_point;
	m_firstTetrahedron = vertex.m_firstTetrahedron;
    m_check            = vertex.m_check;
	m_property         = vertex.m_property;
}



T3DVertex::~T3DVertex()
{
}




rg_Point3D T3DVertex::getPoint() const
{
	return m_point;
}



T3DTetrahedron* T3DVertex::getFirstTetrahedron() const
{
	return m_firstTetrahedron;
}



void* T3DVertex::getProperty() const
{
	return m_property;
}



rg_FLAG T3DVertex::isChecked() const
{
    if ( m_check >= ON_PROCESS )
        return ON_PROCESS;
    else
        return m_check;
}


void T3DVertex::setPoint(const rg_Point3D& point)
{
	m_point = point;
}



void T3DVertex::setFirstTetrahedron(T3DTetrahedron* firstTetrahedron)
{
	m_firstTetrahedron = firstTetrahedron;
}



void T3DVertex::setCheck(const rg_FLAG& check)
{
    if ( m_check >= ON_PROCESS )  {
        if ( check == NOT_VISITED )
            m_check = check;
        else if ( check == ON_PROCESS )
            m_check++;
        else // if ( check == COMPLETED ) 
            m_check--;
    }
    else {
        m_check = check;
    }
}



void T3DVertex::setProperty(void* property)
{
	m_property = property;
}



void T3DVertex::setVertex(const rg_Point3D& point, T3DTetrahedron* firstTetrahedron, void* property)
{
	m_point            = point;
	m_firstTetrahedron = firstTetrahedron;
	m_property         = property;
}




T3DVertex& T3DVertex::operator =(const T3DVertex& vertex)
{
	if ( this == &vertex )
		return *this;

	m_ID               = vertex.m_ID;
	m_point            = vertex.m_point;
	m_firstTetrahedron = vertex.m_firstTetrahedron;
    m_check            = vertex.m_check;
	m_property         = vertex.m_property;

	return *this;
}



