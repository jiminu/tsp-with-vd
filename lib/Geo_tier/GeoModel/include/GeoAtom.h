#ifndef _GEOATOM_H
#define _GEOATOM_H

#include "Sphere.h"

class GeoAtom : public Sphere
{
private:
    rg_INT      m_ID;
	rg_FLOAT    m_RGBColor[3];

public:	
    GeoAtom();
    GeoAtom( const rg_INT& ID );
    GeoAtom( const rg_INT& ID, 
              const rg_REAL& xCoord,  const rg_REAL& yCoord,    const rg_REAL& zCoord,   const rg_REAL& radius,
              const rg_FLOAT& redVal, const rg_FLOAT& greenVal, const rg_FLOAT& blueVal);
    GeoAtom( const rg_INT& ID, const rg_Point3D& center, const rg_REAL& radius,
              const rg_FLOAT& redVal, const rg_FLOAT& greenVal, const rg_FLOAT& blueVal);
    GeoAtom( const rg_INT& ID, const Sphere& sphere,
              const rg_FLOAT& redVal, const rg_FLOAT& greenVal, const rg_FLOAT& blueVal);
    GeoAtom(const GeoAtom& aGeoAtom);
    ~GeoAtom();


    // Get functions    
    rg_INT      getID() const;
    const rg_FLOAT*   getRGBColor() const;
    rg_FLOAT    getRGBColor( const rg_INT& colorID ) const;

    Sphere      geometry() const;

    // Set functions
    void        setID( const rg_INT& ID );
    void        setRGBColor( rg_FLOAT* aRGBColor );
    void        setRGBColor( const rg_FLOAT& redVal, const rg_FLOAT& greenVal, const rg_FLOAT& blueVal );
    void        setRGBColor( const rg_INT& i_color, const rg_FLOAT& colorVal );


    //Operator overloading
    GeoAtom& operator=( const GeoAtom& geoAtom );

};

#endif
