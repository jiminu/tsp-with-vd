#include "GeoAtom.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GeoAtom::GeoAtom() 
: Sphere()
{   
    m_ID          = -1;
    m_RGBColor[0] = 0.0;
    m_RGBColor[1] = 0.0;
    m_RGBColor[2] = 0.0;
}



GeoAtom::GeoAtom(const rg_INT& ID) 
: Sphere()

{
    m_ID          = ID;
    m_RGBColor[0] = 0.0;
    m_RGBColor[1] = 0.0;
    m_RGBColor[2] = 0.0;
}



GeoAtom::GeoAtom(const rg_INT& ID, 
                    const rg_REAL& xCoord, const rg_REAL& yCoord, const rg_REAL& zCoord, 
                    const rg_REAL& radius,
                    const rg_FLOAT& redVal, const rg_FLOAT& greenVal, const rg_FLOAT& blueVal) 
: Sphere(xCoord, yCoord, zCoord, radius), m_ID(ID)
{
    m_RGBColor[0] = redVal;
    m_RGBColor[1] = greenVal;
    m_RGBColor[2] = blueVal;
}



GeoAtom::GeoAtom(const rg_INT& ID, const rg_Point3D& center, const rg_REAL& radius,
                    const rg_FLOAT& redVal, const rg_FLOAT& greenVal, const rg_FLOAT& blueVal) 
: Sphere(center, radius), m_ID(ID)
{
    m_RGBColor[0] = redVal;
    m_RGBColor[1] = greenVal;
    m_RGBColor[2] = blueVal;
}



GeoAtom::GeoAtom(const rg_INT& ID, const Sphere& sphere,
                    const rg_FLOAT& redVal, const rg_FLOAT& greenVal, const rg_FLOAT& blueVal) 
: Sphere(sphere), m_ID(ID)

{
    m_RGBColor[0] = redVal;
    m_RGBColor[1] = greenVal;
    m_RGBColor[2] = blueVal;
}



GeoAtom::GeoAtom( const GeoAtom& geoAtom )
: Sphere( geoAtom ), m_ID(geoAtom.m_ID)
{
    m_RGBColor[0] = geoAtom.m_RGBColor[0];
    m_RGBColor[1] = geoAtom.m_RGBColor[1];
    m_RGBColor[2] = geoAtom.m_RGBColor[2];
}



GeoAtom::~GeoAtom()
{
}



rg_INT GeoAtom::getID() const
{
    return m_ID;    
}



const rg_FLOAT* GeoAtom::getRGBColor() const
{
    return m_RGBColor;
}



rg_FLOAT GeoAtom::getRGBColor(const rg_INT& colorID) const
{
    return m_RGBColor[colorID];
}


Sphere       GeoAtom::geometry() const
{
    return Sphere(m_center, m_radius);
}



void GeoAtom::setID(const rg_INT& ID)
{
    m_ID = ID;
}



void GeoAtom::setRGBColor(rg_FLOAT* RGBColor)
{
    m_RGBColor[0] = RGBColor[0];
    m_RGBColor[1] = RGBColor[1];
    m_RGBColor[2] = RGBColor[2];
}



void GeoAtom::setRGBColor(const rg_FLOAT& redVal, const rg_FLOAT& greenVal, const rg_FLOAT& blueVal)
{
    m_RGBColor[0] = redVal;
    m_RGBColor[1] = greenVal;
    m_RGBColor[2] = blueVal;
}



void GeoAtom::setRGBColor(const rg_INT& i_color, const rg_FLOAT& colorVal)
{
    m_RGBColor[i_color] = colorVal;
}




GeoAtom& GeoAtom::operator=( const GeoAtom& geoAtom )
{
    if( this == &geoAtom )
		return *this;
	
    Sphere::operator=(geoAtom);

    m_ID          = geoAtom.m_ID;

    m_RGBColor[0] = geoAtom.m_RGBColor[0];
    m_RGBColor[1] = geoAtom.m_RGBColor[1];
    m_RGBColor[2] = geoAtom.m_RGBColor[2];

    
	return *this; 
}

