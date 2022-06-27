#include "GeoModel.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GeoModel::GeoModel()
{

}



GeoModel::GeoModel( const GeoModel& geoModel )
{
    m_ID        = geoModel.m_ID;
    m_geoAtoms  = geoModel.m_geoAtoms;
    m_fileName  = geoModel.m_fileName;
}



GeoModel::~GeoModel()
{

}



rg_INT GeoModel::getGeoModelID() const
{
    return m_ID;
}



list< GeoAtom >* GeoModel::getGeoAtoms()
{
    return &m_geoAtoms;
}


const list< GeoAtom >&    GeoModel::geoAtoms() const
{
    return m_geoAtoms;
}



string GeoModel::getFileName() const
{
    return m_fileName;
}



void GeoModel::setGeoAtoms( const list< GeoAtom >& geoAtoms )
{
    m_geoAtoms = geoAtoms;
}



GeoAtom* GeoModel::addGeoAtom( const GeoAtom& geoAtom )
{
    m_geoAtoms.push_back(geoAtom);

    list< GeoAtom >::iterator it_geoAtom = m_geoAtoms.end();

    it_geoAtom--;

    return &(*it_geoAtom);
}



void GeoModel::setFileName( const string& fileName )
{
    m_fileName = fileName;
}



GeoModel& GeoModel::operator=( const GeoModel& geoModel )
{
    if ( this == &geoModel )
        return *this;
    
    m_geoAtoms.clear();

    m_ID        = geoModel.m_ID;
    m_geoAtoms  = geoModel.m_geoAtoms;
    m_fileName  = geoModel.m_fileName;

    return *this;
}

void GeoModel::clear()
{
	m_ID = -1;
	m_fileName = "";
	m_geoAtoms.clear();
}

