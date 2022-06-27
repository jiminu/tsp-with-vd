#include "DiskGenerator2D.h"
#include "VFace2D.h"
using namespace V::GeometryTier;


DiskGenerator2D::DiskGenerator2D()
    : Generator2D()
{
}



DiskGenerator2D::DiskGenerator2D(const int& ID)
    : Generator2D(ID)
{
}



DiskGenerator2D::DiskGenerator2D(const int& ID, const rg_Circle2D& disk)
    : Generator2D(ID, disk)
{
}



DiskGenerator2D::DiskGenerator2D(const int& ID, const rg_Circle2D& disk, VFace2D* const face)
    : Generator2D(ID, disk, face)
{
}



DiskGenerator2D::DiskGenerator2D(const int& ID, const rg_Circle2D& disk, void* const userData)
    : Generator2D(ID, disk, userData)
{
}



DiskGenerator2D::DiskGenerator2D(const DiskGenerator2D& generator)
    : Generator2D(generator)
{
}



DiskGenerator2D::~DiskGenerator2D()
{
}



DiskGenerator2D& DiskGenerator2D::operator=(const DiskGenerator2D& generator)
{
	if (this == &generator) {
		return *this;
	}

    m_type          = generator.m_type;
    m_ID            = generator. m_ID;
    m_disk          = generator. m_disk;
    m_outerFace     = generator. m_outerFace;
    m_innerFace     = generator. m_innerFace;
    m_innerVD       = generator. m_innerVD;
    m_userData      = generator. m_userData;
    m_innerGens     = generator. m_innerGens;

	return *this;
}




bool DiskGenerator2D::operator==(const DiskGenerator2D& generator) const
{
	if (this == &generator) {
		return true;
	}

    if(     m_type          == generator.m_type
        &&  m_ID            == generator.m_ID 
        &&  m_disk          == generator.m_disk 
        &&  m_outerFace     == generator.m_outerFace
        &&  m_innerFace     == generator.m_innerFace
        &&  m_innerVD       == generator.m_innerVD
        &&  m_userData      == generator.m_userData) 
    {
        return true;
    }
    else {
        return false;
    }
}



void DiskGenerator2D::getNeighborGeneratorsInThisVoronoiDiagram(list<DiskGenerator2D*>& neighborGeneratorsList, VoronoiDiagram2DC* const VD)
{
    list<Generator2D*> neighborGens;
    Generator2D::getNeighborGeneratorsInThisVoronoiDiagram(neighborGens, VD);

    for ( list<Generator2D*>::iterator i_gen = neighborGens.begin(); i_gen != neighborGens.end(); ++i_gen ) {
        Generator2D* gen = *i_gen;
        neighborGeneratorsList.push_back((DiskGenerator2D*)gen);
	}
}

void DiskGenerator2D::getNeighborGenerators(list<DiskGenerator2D*>& neighborGeneratorsList, const bool& isThisContainer) const
{
    list<Generator2D*> neighborGens;
    Generator2D::getNeighborGenerators(neighborGens, isThisContainer);

    for (list<Generator2D*>::iterator i_gen = neighborGens.begin(); i_gen != neighborGens.end(); ++i_gen) {
        Generator2D* gen = *i_gen;
        neighborGeneratorsList.push_back((DiskGenerator2D*)gen);
    }
}


