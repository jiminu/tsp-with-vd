#ifndef _GEOMODEL_H
#define _GEOMODEL_H

#include "GeoAtom.h"


#include <list>
#include <string>
using namespace std;

class GeoModel  
{
private:
    rg_INT              m_ID;
    list< GeoAtom >     m_geoAtoms;

    string              m_fileName;

public:
	GeoModel();
    GeoModel(const GeoModel& geoModel);
	~GeoModel();


    // Get functions    
    rg_INT              getGeoModelID() const;
    list< GeoAtom >*    getGeoAtoms();
    string              getFileName() const;

    const list< GeoAtom >&    geoAtoms() const;

    // Set functions
    void                setGeoAtoms(const list< GeoAtom >& geoAtoms);
    GeoAtom*            addGeoAtom(const GeoAtom& geoAtom);
    void                setFileName(const string& fileName);


	// Clear function
	void clear();

    //Operator overloading
    GeoModel& operator =(const GeoModel& geoModel);

};

#endif
