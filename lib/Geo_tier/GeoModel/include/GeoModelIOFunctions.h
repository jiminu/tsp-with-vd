#ifndef _GEOMODEL_IO_FUNCTIONS_H
#define _GEOMODEL_IO_FUNCTIONS_H


#include "rg_Const.h"
#include "rg_dList.h"
#include "GeoModel.h"

#include <fstream>
#include <string>
using namespace std;



namespace GeoModelIOFunctions
{
    
    // MCO IO Functions
	rg_BOOL readFile( const char* fileName, GeoModel& aGeoModel );

};

#endif
