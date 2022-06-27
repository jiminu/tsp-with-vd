#include "GeoModelIOFunctions.h"

rg_BOOL GeoModelIOFunctions::readFile( const char* fileName, GeoModel& aGeoModel )
{
	rg_FLAG isFileRead = rg_TRUE;

    char        lineOfData[200];
    const char* seps    = " \t\n\r";
    char*       token   = NULL;

    ifstream fin;
    fin.open(fileName);           

    string startChar;
    while ( !fin.eof() ) {
        fin.getline(lineOfData, 200);
        startChar = lineOfData[0];
        if( startChar != "%" && startChar != "#" )
            break;
    }

    token = strtok( lineOfData, seps );
    rg_INT numOfGeoAtoms = atoi( token );
    
    rg_INT  id       = 0;
    rg_REAL xCoord   = 0.0;
    rg_REAL yCoord   = 0.0;
    rg_REAL zCoord   = 0.0;
    rg_REAL radius   = 0.0;
    rg_REAL redVal   = 0.0;
    rg_REAL greenVal = 0.0;
    rg_REAL blueVal  = 0.0;

//     for( rg_INT i_particle=0; i_particle<numOfGeoAtoms; i_particle++ ) {

    rg_INT countGeoAtom = 0;
    while ( !fin.eof() && countGeoAtom < numOfGeoAtoms) {
        fin.getline(lineOfData, 200);

        startChar = lineOfData[0];
        if( startChar == "%" || startChar == "#" ) {
            continue;
        }


        // check if record line does not exists!!

//         fin.getline(lineOfData, 200);

        if( !strcmp( lineOfData, "" ) ) {
            isFileRead = rg_FALSE;
            aGeoModel.clear();
            break;
        }
        
        // ID
        token    = strtok( lineOfData, seps );
        id       = atoi( token );

        // Coordinates
        xCoord   = atof( strtok( NULL, seps ) );
        yCoord   = atof( strtok( NULL, seps ) );
        zCoord   = atof( strtok( NULL, seps ) );

        // Radius
        radius   = atof( strtok( NULL, seps ) );

        // RGB colors
        token    = strtok( NULL, seps );

        if( token == NULL ) {
            redVal   = 0.7;
            greenVal = 0.7;
            blueVal  = 0.7;            
        }
        else {
            redVal   = atof( token );
            greenVal = atof( strtok( NULL, seps ) );
            blueVal  = atof( strtok( NULL, seps ) );
        }

        aGeoModel.addGeoAtom( GeoAtom(id, xCoord, yCoord, zCoord, radius, redVal, greenVal, blueVal) ); 
        countGeoAtom++;
    }

	aGeoModel.setFileName(fileName);

    fin.close();

    if ( countGeoAtom == numOfGeoAtoms ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}

