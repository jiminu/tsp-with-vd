#include "rg_QTFReader.h"
#include "ConstForQuasiTriangulation.h"
using namespace V::GeometryTier;

#include "StringFunctions.h"
#include <stdlib.h>
#include "rg_dList.h"
#include "Sphere.h"
#include <time.h>

#include <strstream>
using namespace std;


QTFReader::QTFReader()
: m_QT(rg_NULL), m_vertices(rg_NULL), m_cells(rg_NULL), m_useVDVerticesInQTFile(rg_FALSE)
{
}



QTFReader::QTFReader(QuasiTriangulation* QT)
: m_QT(QT), m_vertices(rg_NULL), m_cells(rg_NULL), m_useVDVerticesInQTFile(rg_FALSE)
{
}



QTFReader::QTFReader(const QTFReader& reader)
: m_QT(reader.m_QT), m_vertices(reader.m_vertices), m_cells(reader.m_cells), 
  m_useVDVerticesInQTFile(reader.m_useVDVerticesInQTFile)
{
}



QTFReader::~QTFReader()
{
    m_QT = rg_NULL;
    if ( m_vertices != rg_NULL ) {
        delete [] m_vertices;
    }

    if ( m_cells != rg_NULL ) {
        delete [] m_cells;
    }
}




QuasiTriangulation* QTFReader::getQuasiTriangulation() const
{
    return m_QT;
}



void QTFReader::setQuasiTriangulation(QuasiTriangulation* QT)
{
    m_QT = QT;
}




QTFReader& QTFReader::operator=(const QTFReader& reader)
{
    if ( this == &reader) {
        return *this;
    }

    if ( m_vertices != rg_NULL ) {
        delete [] m_vertices;
    }

    if ( m_cells != rg_NULL ) {
        delete [] m_cells;
    }

    m_QT = reader.m_QT;
    m_vertices = reader.m_vertices;
    m_cells    = reader.m_cells;
    m_useVDVerticesInQTFile = reader.m_useVDVerticesInQTFile;


    return *this;
}


rg_BOOL QTFReader::read(const Molecule* const& molecule, const string& qtfFilename, QuasiTriangulation* QT, const rg_BOOL& useVDVertices)
{
    m_QT = QT;


    rg_BOOL isValidQTF = rg_TRUE;

    //  read QTF header and validate QTF
    string pdbTimeStampOfMolecule = molecule->getTimeStamp();
    rg_INT pdbFileSizeOfMolecule  = molecule->getFileSize();

    string pdbTimeStampInQTF      = parsePDBTimeStamp( qtfFilename );
    rg_INT pdbFileSizeInQTF       = parsePDBFileSize(  qtfFilename );

    if ( (pdbTimeStampOfMolecule != pdbTimeStampInQTF) || (pdbFileSizeOfMolecule != pdbFileSizeInQTF) ) {
        return rg_FALSE;
    }


    //  read QTF body
    read(qtfFilename, molecule->getModelSerialFromInputFile(), useVDVertices);


    //  synchronize atoms between molecule and QTF
    isValidQTF = synchronizeAtomsBetweenQTAndMolecule( molecule );


    return isValidQTF;
}



string QTFReader::parsePDBTimeStamp(const string& qtfFilename)
{
    ifstream fin;
    fin.open( qtfFilename.c_str() );

    string timeStamp;

    string strRecord  = "";
    while ( getline( fin, strRecord ) ) {
        rg_INT position = 0;
        string recType  = StringFunctions::subString( strRecord, position );

        if( recType == QTF_RECORD_PDBTIMESTAMP ) {
            timeStamp = StringFunctions::subString( strRecord, position );
            break;
        }
    }

    fin.close();

    return timeStamp;
}



rg_INT QTFReader::parsePDBFileSize(const string& qtfFilename)
{
    ifstream fin;
    fin.open( qtfFilename.c_str() );

    rg_INT pdbFileSize = 0;

    string strRecord  = "";
    while ( getline( fin, strRecord ) ) {
        rg_INT position = 0;
        string recType  = StringFunctions::subString( strRecord, position );

        if( recType == QTF_RECORD_PDBFILESIZE ) {
            pdbFileSize = atoi( StringFunctions::subString( strRecord, position ).c_str() );
            break;
        }
    }

    fin.close();

    return pdbFileSize;
}



rg_BOOL QTFReader::read(const string& qtfFilename, const rg_INT& modelID, const rg_BOOL& useVDVertices)
{
    m_useVDVerticesInQTFile = useVDVertices;

    list<string> recordsForQTVertices;
    list<string> recordsForQTCells;
    list<string> recordsForQTGates;
    list<string> recordsForVDVertices;

    extractSingleModelFromQTFile(qtfFilename, modelID, recordsForQTVertices, recordsForQTCells, recordsForQTGates, recordsForVDVertices);

    parseVerticesOfQTFile( recordsForQTVertices );
    parseCellsOfQTFile(    recordsForQTCells );

    setTopologyBetweenVerticesAndCells();

    parseGatesOfQTFile(    recordsForQTGates );

    if ( m_useVDVerticesInQTFile ) {
        parseVDVerticesOfQTFile(  recordsForVDVertices );
    }
    else {
        computeEmptyTangentSpheresOfCells();
    }
    
    return rg_TRUE;
}

rg_BOOL QTFReader::read2(const string& qtfFilename, const rg_INT& modelID, const rg_BOOL& useVDVertices)
{
    
    m_useVDVerticesInQTFile = useVDVertices;


    list<string> recordsOfSelectedModel;
    extractSingleModelFromQTFile(qtfFilename, modelID, recordsOfSelectedModel );


    list<string*> recordsForCells;
    list<string*> recordsForGates;
    parseVerticesOfQTFile(recordsOfSelectedModel, recordsForCells );
    parseCellsOfQTFile(   recordsForCells,        recordsForGates );


    setTopologyBetweenVerticesAndCells();


    list<string*> recordsForVDVertices;
    parseGatesOfQTFile(   recordsForGates,        recordsForVDVertices);

    if ( m_useVDVerticesInQTFile ) {
        parseVDVerticesOfQTFile( recordsForVDVertices );
    }
    else {
        computeEmptyTangentSpheresOfCells();
    }


    return rg_TRUE;
    


    /*
    // for analysis of QTFReader
	clock_t start, finish;   
    timeLoading = 0.0;
    timeParsing = 0.0;
    timeSettingQT = 0.0;
    timeParsingVVtx = 0.0;
    timeComputingVVtx = 0.0;

    m_useVDVerticesInQTFile = useVDVertices;


    list<string> recordsOfSelectedModel;

	start = clock();
        extractSingleModelFromQTFile(qtfFilename, modelID, recordsOfSelectedModel );
	finish = clock();
	timeLoading = (double)(finish - start) / CLOCKS_PER_SEC;


    list<string*> recordsForCells;
    list<string*> recordsForGates;

	start = clock();
        parseVerticesOfQTFile(recordsOfSelectedModel, recordsForCells );
        parseCellsOfQTFile(   recordsForCells,        recordsForGates );
	finish = clock();
	timeParsing = (double)(finish - start) / CLOCKS_PER_SEC;


	start = clock();
        setTopologyBetweenVerticesAndCells();
	finish = clock();
	timeSettingQT = (double)(finish - start) / CLOCKS_PER_SEC;


    list<string*> recordsForVDVertices;
    parseGatesOfQTFile(   recordsForGates,        recordsForVDVertices);

    if ( m_useVDVerticesInQTFile ) {
	    start = clock();
        parseVDVerticesOfQTFile( recordsForVDVertices );
	    finish = clock();
	    timeParsingVVtx = (double)(finish - start) / CLOCKS_PER_SEC;
    }
    else {
	    start = clock();
        computeEmptyTangentSpheresOfCells();
	    finish = clock();
	    timeComputingVVtx = (double)(finish - start) / CLOCKS_PER_SEC;
    }


    return rg_TRUE;
    */
}



void QTFReader::extractSingleModelFromQTFile(const string& qtfFilename, const rg_INT& modelID, list<string>& recordsOfSelectedModel )
{
    ifstream fin;
    fin.open( qtfFilename.c_str() );

    fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE
    
    rg_BOOL isRecordInSelectedModel = rg_FALSE;

    string strRecord  = "";
    while ( getline( fin, strRecord ) ) {
        rg_INT position = 0;
        string recType  = StringFunctions::subString( strRecord, position );

        if ( !isRecordInSelectedModel ) {
            if( recType == QTF_RECORD_MODEL ) {
                rg_INT  currModelID = atoi( StringFunctions::subString( strRecord, position ).c_str() );

                if ( modelID == currModelID ) {
                    isRecordInSelectedModel = rg_TRUE;
                }
            }
        }   
        else {
            if( recType == QTF_RECORD_ENDMODEL ) {
                break;
            }

            recordsOfSelectedModel.push_back( StringFunctions::strTrimRight( strRecord ) );
        }
    }

    fin.close();
}



void QTFReader::extractSingleModelFromQTFile( const string& qtfFilename, const rg_INT& modelID, 
                                              list<string>& recordsForQTVertices,
                                              list<string>& recordsForQTCells,
                                              list<string>& recordsForQTGates,
                                              list<string>& recordsForVDVertices )
{
    ifstream fin;
    fin.open( qtfFilename.c_str() );

    fin.seekg(0, ios::beg); // GOTO BEGINING OF FILE
    
    rg_BOOL isRecordInSelectedModel = rg_FALSE;

    string strRecord  = "";
    while ( getline( fin, strRecord ) ) {
        rg_INT position = 0;
        string recType  = StringFunctions::subString( strRecord, position );

        if ( !isRecordInSelectedModel ) {
            if( recType == QTF_RECORD_MODEL ) {
                rg_INT  currModelID = atoi( StringFunctions::subString( strRecord, position ).c_str() );

                if ( modelID == currModelID ) {
                    isRecordInSelectedModel = rg_TRUE;
                }
            }
            continue;
        }   

        if( recType == QTF_RECORD_VERTEX ) {
            recordsForQTVertices.push_back( StringFunctions::strTrimRight( strRecord ) );
        }
        else if( recType == QTF_RECORD_CELL ) {
            recordsForQTCells.push_back( StringFunctions::strTrimRight( strRecord ) );
        }
        else if( recType == QTF_RECORD_GATE ) {
            recordsForQTGates.push_back( StringFunctions::strTrimRight( strRecord ) );
        }
        else if( recType == QTF_RECORD_VDVTX ) {
            recordsForVDVertices.push_back( StringFunctions::strTrimRight( strRecord ) );
        }
        else if( recType == QTF_RECORD_ENDMODEL ) {
            break;
        }

    }

    fin.close();
}



void QTFReader::parseVerticesOfQTFile(list<string>& recordsInModel, list<string*>& recordsForCells )
{
    char  seps[]   = " ,\t\n";

    m_QT->addVertex( QTVertex(0, rg_NULL) );


    list<string>::iterator i_recLines;
    for( i_recLines = recordsInModel.begin(); i_recLines != recordsInModel.end(); i_recLines++ ) {
        string currRecLine = (*i_recLines);

        //  using strtok() - 2012.02.27 by Y.Cho
        const char* record = currRecLine.c_str();

        string recType   = strtok( (char*)record, seps );
        if( recType != QTF_RECORD_VERTEX ) {
            recordsForCells.push_back( &(*i_recLines) );
            continue;
        }

        rg_INT  vertexID    = atoi( strtok( NULL, seps ) );
        rg_INT  firstCellID = atoi( strtok( NULL, seps ) );
        rg_REAL x           = atof( strtok( NULL, seps ) );
        rg_REAL y           = atof( strtok( NULL, seps ) );
        rg_REAL z           = atof( strtok( NULL, seps ) );
        rg_REAL radius      = atof( strtok( NULL, seps ) );
        rg_INT  IDFromInput = atoi( strtok( NULL, seps ) );

        Ball*     ball   = m_QT->addBall( Ball( vertexID, Sphere(x, y, z, radius), rg_NULL, IDFromInput) );
        QTVertex* vertex = m_QT->addVertex( QTVertex(vertexID, ball) );

        QTIncexVertex indexVertex;
        indexVertex.qtVertex  = vertex;
        indexVertex.firstCell = firstCellID;
        m_indexVertices.add( indexVertex );


        /*
        //  original codes for QTFReader
        rg_INT position = 0;
        string recType  = StringFunctions::subString( currRecLine, position );

        if( recType != QTF_RECORD_VERTEX ) {
            recordsForCells.push_back( &(*i_recLines) );
            continue;
        }

        rg_INT  vertexID    = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        rg_INT  firstCellID = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        rg_REAL x           = atof( StringFunctions::subString( currRecLine, position ).c_str() );
        rg_REAL y           = atof( StringFunctions::subString( currRecLine, position ).c_str() );
        rg_REAL z           = atof( StringFunctions::subString( currRecLine, position ).c_str() );
        rg_REAL radius      = atof( StringFunctions::subString( currRecLine, position ).c_str() );
        rg_INT  IDFromInput = atoi( StringFunctions::subString( currRecLine, position ).c_str() );

        Ball*     ball   = m_QT->addBall( Ball( vertexID, Sphere(x, y, z, radius), rg_NULL, IDFromInput) );
        QTVertex* vertex = m_QT->addVertex( QTVertex(vertexID, ball) );

        QTIncexVertex indexVertex;
        indexVertex.qtVertex  = vertex;
        indexVertex.firstCell = firstCellID;
        m_indexVertices.add( indexVertex );
        */


        // 아래의 코드(strstream)가 위의 것(StringFunctions)보다 4배정도 더 느림.
        //istrstream strStream( currRecLine.c_str() );

        //string recType;
        //strStream >> recType;
        //if( recType != QTF_RECORD_VERTEX ) {
        //    recordsForCells.push_back( &(*i_recLines) );
        //    continue;
        //}

        //rg_INT  vertexID    = 0;
        //rg_INT  firstCellID = 0;
        //rg_REAL x           = 0.0;
        //rg_REAL y           = 0.0;
        //rg_REAL z           = 0.0;
        //rg_REAL radius      = 0.0;
        //rg_INT  IDFromInput = 0;

        //strStream >> vertexID >> firstCellID >> x >> y >> z >> radius >> IDFromInput; 

        //Ball*     ball   = m_QT->addBall( Ball( vertexID, Sphere(x, y, z, radius), rg_NULL, IDFromInput) );
        //QTVertex* vertex = m_QT->addVertex( QTVertex(vertexID, ball) );

        //QTIncexVertex indexVertex;
        //indexVertex.qtVertex  = vertex;
        //indexVertex.firstCell = firstCellID;
        //m_indexVertices.add( indexVertex );

        /*
        //  using StringFunctions::strToken()
        rg_INT position = 0;
        string recType  = StringFunctions::strToken( currRecLine, position );

        if( recType != QTF_RECORD_VERTEX ) {
            recordsForCells.push_back( &(*i_recLines) );
            continue;
        }

        rg_INT  vertexID    = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_INT  firstCellID = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_REAL x           = atof( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_REAL y           = atof( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_REAL z           = atof( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_REAL radius      = atof( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_INT  IDFromInput = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );

        Ball*     ball   = m_QT->addBall( Ball( vertexID, Sphere(x, y, z, radius), rg_NULL, IDFromInput) );
        QTVertex* vertex = m_QT->addVertex( QTVertex(vertexID, ball) );

        QTIncexVertex indexVertex;
        indexVertex.qtVertex  = vertex;
        indexVertex.firstCell = firstCellID;
        m_indexVertices.add( indexVertex );
        */
    }
}

void QTFReader::parseVerticesOfQTFile(const list<string>& recordsForQTVertices )
{
    char  seps[]   = " ,\t\n";

    m_QT->addVertex( QTVertex(0, rg_NULL) );


    list<string>::const_iterator i_recLines;
    for( i_recLines = recordsForQTVertices.begin(); i_recLines != recordsForQTVertices.end(); i_recLines++ ) {
        string currRecLine = (*i_recLines);

        //  using strtok() - 2012.02.27 by Y.Cho
        const char* record = currRecLine.c_str();

        string  recType     = strtok( (char*)record, seps );

        rg_INT  vertexID    = atoi( strtok( NULL, seps ) );
        rg_INT  firstCellID = atoi( strtok( NULL, seps ) );
        rg_REAL x           = atof( strtok( NULL, seps ) );
        rg_REAL y           = atof( strtok( NULL, seps ) );
        rg_REAL z           = atof( strtok( NULL, seps ) );
        rg_REAL radius      = atof( strtok( NULL, seps ) );
        rg_INT  IDFromInput = atoi( strtok( NULL, seps ) );

        Ball*     ball   = m_QT->addBall( Ball( vertexID, Sphere(x, y, z, radius), rg_NULL, IDFromInput) );
        QTVertex* vertex = m_QT->addVertex( QTVertex(vertexID, ball) );

        QTIncexVertex indexVertex;
        indexVertex.qtVertex  = vertex;
        indexVertex.firstCell = firstCellID;
        m_indexVertices.add( indexVertex );
    }
}



void QTFReader::parseCellsOfQTFile(   const list<string>& recordsForQTCells )
{
    char  seps[]   = " ,\t\n";

    list<string>::const_iterator i_recLines = recordsForQTCells.begin();
    for( ; i_recLines != recordsForQTCells.end(); i_recLines++ ) {
        string currRecLine = (*i_recLines);

        //  using strtok() - 2012.02.27 by Y.Cho
        const char* record = currRecLine.c_str();

        string recType   = strtok( (char*)record, seps );

        rg_INT  cellID   = atoi( strtok( NULL, seps ) );

        QTIndexCell indexCell;
        indexCell.vertexIndex[0]   = atoi( strtok( NULL, seps ) );
        indexCell.vertexIndex[1]   = atoi( strtok( NULL, seps ) );
        indexCell.vertexIndex[2]   = atoi( strtok( NULL, seps ) );
        indexCell.vertexIndex[3]   = atoi( strtok( NULL, seps ) );
        indexCell.neighborIndex[0] = atoi( strtok( NULL, seps ) );
        indexCell.neighborIndex[1] = atoi( strtok( NULL, seps ) );
        indexCell.neighborIndex[2] = atoi( strtok( NULL, seps ) );
        indexCell.neighborIndex[3] = atoi( strtok( NULL, seps ) );

        indexCell.qtCell = m_QT->addTetrahedron( QTTetrahedron( cellID-1 ) );
        m_indexCells.add( indexCell );
    }
}



void QTFReader::parseGatesOfQTFile(   const list<string>& recordsForQTGates )
{
    char* strToken;
    char  seps[]   = " ,\t\n";

    list<string>::const_iterator i_recLines;
    for( i_recLines = recordsForQTGates.begin(); i_recLines != recordsForQTGates.end(); i_recLines++ ) {
        string currRecLine = (*i_recLines);

        const char* record = currRecLine.c_str();

        string recType   = strtok( (char*)record, seps );

        rg_INT gateID    = atoi( strtok( NULL, seps ) );
        rg_INT vertexID1 = atoi( strtok( NULL, seps ) );
        rg_INT vertexID2 = atoi( strtok( NULL, seps ) );
        rg_INT bigCellID = atoi( strtok( NULL, seps ) );

        QTGate* currGate = m_QT->addGate( QTGate(m_vertices[vertexID1], m_vertices[vertexID2]) );

        currGate->setBigTetrahedron( m_cells[bigCellID] );
    

        strToken = strtok( NULL, seps );
        while ( strToken != rg_NULL ) {
            rg_INT smallWorldID = atoi( strToken );
            currGate->addSmallTetrahedron( m_cells[smallWorldID] );

            strToken = strtok( NULL, seps );
        }
    }
}


void QTFReader::parseVDVerticesOfQTFile(const list<string>& recordsForVDVertices)
{
    char  seps[]   = " ,\t\n";

    list<string>::const_iterator i_recLines;
    for( i_recLines = recordsForVDVertices.begin(); i_recLines != recordsForVDVertices.end(); i_recLines++ ) {
        string currRecLine = (*i_recLines);

        //  using strtok() - 2012.02.27 by Y.Cho
        const char* record = currRecLine.c_str();

        string recType   = strtok( (char*)record, seps );

        rg_INT  cellID = atoi( strtok( NULL, seps ) );
        rg_REAL x      = atof( strtok( NULL, seps ) );
        rg_REAL y      = atof( strtok( NULL, seps ) );
        rg_REAL z      = atof( strtok( NULL, seps ) );

        Sphere  atom  = m_cells[cellID]->getVertex(0)->getBall();

        rg_Point3D center(x, y, z);
        rg_REAL    radius = center.distance( atom.getCenter() ) - atom.getRadius();

        m_cells[cellID]->setEmptyTangentSphere( Sphere(x, y, z, radius) );
    }
}




void QTFReader::parseCellsOfQTFile(   list<string*>& recordsForCells, list<string*>& recordsForGates )
{
    char  seps[]   = " ,\t\n";

    list<string*>::iterator i_recLines = recordsForCells.begin();
    for( ; i_recLines != recordsForCells.end(); i_recLines++ ) {
        string currRecLine = *(*i_recLines);

        //  using strtok() - 2012.02.27 by Y.Cho
        const char* record = currRecLine.c_str();

        string recType   = strtok( (char*)record, seps );
        if( recType != QTF_RECORD_CELL ) {
            recordsForGates.push_back((*i_recLines));
            continue;
        }

        rg_INT  cellID = atoi( strtok( NULL, seps ) );

        QTIndexCell indexCell;
        indexCell.vertexIndex[0]   = atoi( strtok( NULL, seps ) );
        indexCell.vertexIndex[1]   = atoi( strtok( NULL, seps ) );
        indexCell.vertexIndex[2]   = atoi( strtok( NULL, seps ) );
        indexCell.vertexIndex[3]   = atoi( strtok( NULL, seps ) );
        indexCell.neighborIndex[0] = atoi( strtok( NULL, seps ) );
        indexCell.neighborIndex[1] = atoi( strtok( NULL, seps ) );
        indexCell.neighborIndex[2] = atoi( strtok( NULL, seps ) );
        indexCell.neighborIndex[3] = atoi( strtok( NULL, seps ) );

        indexCell.qtCell = m_QT->addTetrahedron( QTTetrahedron( cellID-1 ) );
        m_indexCells.add( indexCell );
        
        
        /*
        //  original codes for QTFReader
        rg_INT position = 0;
        string recType  = StringFunctions::subString( currRecLine, position );

        if( recType != QTF_RECORD_CELL ) {
            recordsForGates.push_back((*i_recLines));
            continue;
        }

        rg_INT  cellID = atoi( StringFunctions::subString( currRecLine, position ).c_str() );

        QTIndexCell indexCell;
        indexCell.vertexIndex[0]   = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        indexCell.vertexIndex[1]   = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        indexCell.vertexIndex[2]   = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        indexCell.vertexIndex[3]   = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        indexCell.neighborIndex[0] = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        indexCell.neighborIndex[1] = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        indexCell.neighborIndex[2] = atoi( StringFunctions::subString( currRecLine, position ).c_str() );
        indexCell.neighborIndex[3] = atoi( StringFunctions::subString( currRecLine, position ).c_str() );

        indexCell.qtCell = m_QT->addTetrahedron( QTTetrahedron( cellID-1 ) );
        m_indexCells.add( indexCell );
        */


        // 아래의 코드(strstream)가 위의 것(StringFunctions)보다 4배정도 더 느림.
        //istrstream strStream( currRecLine.c_str() );

        //string recType;
        //strStream >> recType;
        //if( recType != QTF_RECORD_CELL ) {
        //    recordsForGates.push_back((*i_recLines));
        //    continue;
        //}

        //rg_INT      cellID    = 0;
        //QTIndexCell indexCell;

        //strStream >> cellID 
        //    >> indexCell.vertexIndex[0]   >> indexCell.vertexIndex[1]   >> indexCell.vertexIndex[2]   >> indexCell.vertexIndex[3]
        //    >> indexCell.neighborIndex[0] >> indexCell.neighborIndex[1] >> indexCell.neighborIndex[2] >> indexCell.neighborIndex[3];

        //indexCell.qtCell = m_QT->addTetrahedron( QTTetrahedron( cellID-1 ) );
        //m_indexCells.add( indexCell );

        /*
        //  using StringFunctions::strToken()
        rg_INT position = 0;
        string recType  = StringFunctions::strToken( currRecLine, position );

        if( recType != QTF_RECORD_CELL ) {
            recordsForGates.push_back((*i_recLines));
            continue;
        }

        rg_INT  cellID = atoi( StringFunctions::subString( currRecLine, position ).c_str() );

        QTIndexCell indexCell;
        indexCell.vertexIndex[0]   = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        indexCell.vertexIndex[1]   = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        indexCell.vertexIndex[2]   = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        indexCell.vertexIndex[3]   = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        indexCell.neighborIndex[0] = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        indexCell.neighborIndex[1] = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        indexCell.neighborIndex[2] = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        indexCell.neighborIndex[3] = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );

        indexCell.qtCell = m_QT->addTetrahedron( QTTetrahedron( cellID-1 ) );
        m_indexCells.add( indexCell );
        */

    }
}

    

void QTFReader::setTopologyBetweenVerticesAndCells()
{    
    rg_INT numVertices = m_QT->getNumQTVertices();
    m_vertices = new QTVertex*[ numVertices ];

    rg_dList<QTVertex>* vertexList = m_QT->getVertices();
    vertexList->reset4Loop();
    while ( vertexList->setNext4Loop() ) {
        QTVertex* currVtx = vertexList->getpEntity();

        m_vertices[ currVtx->getID() ] = currVtx;
    }


    rg_INT numCells = m_QT->getNumQTTetrahedron() + 1;
    m_cells = new QTTetrahedron*[ numCells ];

    m_cells[0] = rg_NULL;
    rg_dList<QTTetrahedron>* cellList = m_QT->getTetrahedra();
    cellList->reset4Loop();
    while ( cellList->setNext4Loop() ) {
        QTTetrahedron* currCell = cellList->getpEntity();
        m_cells[ currCell->getID()+1 ] = currCell;
    }


    //  set first q-cell of each q-vertex.
    m_indexVertices.reset4Loop();
    while ( m_indexVertices.setNext4Loop() ) {
        QTIncexVertex currIndexVtx = m_indexVertices.getEntity();

        currIndexVtx.qtVertex->setFirstTetrahedron( m_cells[currIndexVtx.firstCell] );
    }



    //  set four bounding q-vertices and four neighbor q-cells of each q-cell.
    rg_dList<QTTetrahedron*> cellsToContributeToBoundaryOfSqueezedHull;

    m_indexCells.reset4Loop();
    while ( m_indexCells.setNext4Loop() ) {
        QTIndexCell*   currIndexCell = m_indexCells.getpEntity();

        QTTetrahedron* currCell = currIndexCell->qtCell;

        if ( currIndexCell->vertexIndex[3] != -1 ) {
            currCell->setVertices(  m_vertices[ currIndexCell->vertexIndex[0] ], 
                                    m_vertices[ currIndexCell->vertexIndex[1] ], 
                                    m_vertices[ currIndexCell->vertexIndex[2] ], 
                                    m_vertices[ currIndexCell->vertexIndex[3] ] );

            currCell->setNeighbors( m_cells[  currIndexCell->neighborIndex[0] ], 
                                    m_cells[  currIndexCell->neighborIndex[1] ], 
                                    m_cells[  currIndexCell->neighborIndex[2] ], 
                                    m_cells[  currIndexCell->neighborIndex[3] ] );


            // currIndexCell->neighborIndex[i]==0 means 
            // that the corresponding q-face is on the boundary of the squeezed hull.
            if (    currIndexCell->neighborIndex[0] == 0 || currIndexCell->neighborIndex[1] == 0 
                 || currIndexCell->neighborIndex[2] == 0 || currIndexCell->neighborIndex[3] == 0 ) {
                cellsToContributeToBoundaryOfSqueezedHull.add( currCell );
            }
        }
        else {
            //  for an isolated triangle.
            currCell->setVertices(  m_vertices[ currIndexCell->vertexIndex[0] ], 
                                    m_vertices[ currIndexCell->vertexIndex[1] ], 
                                    m_vertices[ currIndexCell->vertexIndex[2] ], 
                                    rg_NULL );
            currCell->setNeighbors( rg_NULL, rg_NULL, rg_NULL, rg_NULL );
        }
    }



    //  set topology of cells incident to virtual vertex(infinite vertex).
    QTVertex* virtualVertex = m_vertices[0];
    rg_dList<QTTetrahedron*> cellsIncidentToVirtualVertex;

    //  make cells incident to virtual vertex.
    cellsToContributeToBoundaryOfSqueezedHull.reset4Loop();
    while ( cellsToContributeToBoundaryOfSqueezedHull.setNext4Loop() ) {
        QTTetrahedron* currCell = cellsToContributeToBoundaryOfSqueezedHull.getEntity();
        
        for ( rg_INT i=0; i<4; i++ ) {
            QTTetrahedron* neighbor = currCell->getNeighbor(i);
            if ( neighbor != rg_NULL ) {
                continue;
            }

            QTFace face( currCell, currCell->getVertex(i) );
            QTVertex* vertex[3] = {rg_NULL, rg_NULL, rg_NULL};
            face.inquireBoundingVerticesInThisWorld( vertex );

            neighbor = m_QT->addTetrahedron( QTTetrahedron( m_QT->getNumQTTetrahedron() ) );
            neighbor->setVertices(virtualVertex, vertex[2], vertex[1], vertex[0] );

            neighbor->setNeighbor(0, currCell);
            currCell->setNeighbor(i, neighbor);

            cellsIncidentToVirtualVertex.add(neighbor);
        }
    }


    //  set topology among cells incident to virtual vertex.
    virtualVertex->setFirstTetrahedron( cellsIncidentToVirtualVertex.getFirstEntity() );

    cellsIncidentToVirtualVertex.reset4Loop();
    while ( cellsIncidentToVirtualVertex.setNext4Loop() ) {
        QTTetrahedron* currCell = cellsIncidentToVirtualVertex.getEntity();

        QTFace faceOnSqueezedHull( currCell, currCell->getVertex(0) );

        rg_INT i=0;
        QTVertex* vertex[6] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL, rg_NULL, rg_NULL};
        for ( i=0; i<4; i++ ) {
            vertex[i] = currCell->getVertex(i);
        }
        vertex[4] = vertex[1];
        vertex[5] = vertex[2];

        for ( i=1; i<4; i++ ) {
            QTTetrahedron* neighbor = currCell->getNeighbor(i);
            if ( neighbor != rg_NULL ) {
                continue;
            }

            QTEdge edge(currCell, vertex[i+1], vertex[i+2]);

            QTTetrahedron* prevCell = rg_NULL;
            QTTetrahedron* nextCell = currCell;
            do {
                prevCell = nextCell;
                nextCell = edge.getCWNextCell( prevCell );
            } while ( nextCell != rg_NULL );

            neighbor = prevCell;
            currCell->setNeighbor(i, neighbor);

            rg_INT posOfCurrCellInNeighbor = -1;
            for ( rg_INT j=1; j<4; j++ ) {
                QTVertex* vetexInNeighbor = neighbor->getVertex(j);
                if ( vetexInNeighbor != vertex[i+1] && vetexInNeighbor != vertex[i+2] ) {
                    posOfCurrCellInNeighbor = j;
                    break;
                }
            }

            neighbor->setNeighbor(posOfCurrCellInNeighbor, currCell);
        }
    }
}



void QTFReader::parseGatesOfQTFile(   list<string*>& recordsForGates, list<string*>& recordsForVDVertices)
{
    list<string*>::iterator i_recLines;
    for( i_recLines = recordsForGates.begin(); i_recLines != recordsForGates.end(); i_recLines++ ) {
        string* currRecLine = (*i_recLines);

        rg_INT position = 0;
        string recType  = StringFunctions::subString( *currRecLine, position );

        if( recType != QTF_RECORD_GATE ) {
            recordsForVDVertices.push_back( currRecLine );
            continue;
        }

        rg_INT gateID    = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_INT vertexID1 = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_INT vertexID2 = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );

        QTGate* currGate = m_QT->addGate( QTGate(m_vertices[vertexID1], m_vertices[vertexID2]) );

        rg_INT bigCellID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        currGate->setBigTetrahedron( m_cells[bigCellID] );
    

        while ( position != -1 ) {
            rg_INT smallWorldID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
            currGate->addSmallTetrahedron( m_cells[smallWorldID] );
        }
    }
}



void QTFReader::parseVDVerticesOfQTFile(list<string*>& recordsForVDVertices)
{
    char  seps[]   = " ,\t\n";

    list<string*>::iterator i_recLines;
    for( i_recLines = recordsForVDVertices.begin(); i_recLines != recordsForVDVertices.end(); i_recLines++ ) {
        string* currRecLine = (*i_recLines);

        //  using strtok() - 2012.02.27 by Y.Cho
        const char* record = currRecLine->c_str();

        string recType   = strtok( (char*)record, seps );
        if( recType != QTF_RECORD_VDVTX ) {
            continue;
        }

        rg_INT  cellID = atoi( strtok( NULL, seps ) );
        rg_REAL x      = atof( strtok( NULL, seps ) );
        rg_REAL y      = atof( strtok( NULL, seps ) );
        rg_REAL z      = atof( strtok( NULL, seps ) );

        Sphere  atom  = m_cells[cellID]->getVertex(0)->getBall();

        rg_Point3D center(x, y, z);
        rg_REAL    radius = center.distance( atom.getCenter() ) - atom.getRadius();

        m_cells[cellID]->setEmptyTangentSphere( Sphere(x, y, z, radius) );
        
        
        /*
        //  original codes for QTFReader

        rg_INT position = 0;
        string recType  = StringFunctions::subString( *currRecLine, position );

        if( recType != QTF_RECORD_VDVTX ) {
            continue;
        }

        rg_INT cellID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_REAL x     = atof( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_REAL y     = atof( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_REAL z     = atof( StringFunctions::subString( *currRecLine, position ).c_str() );

        Sphere  atom  = m_cells[cellID]->getVertex(0)->getBall();

        rg_Point3D center(x, y, z);
        rg_REAL    radius = center.distance( atom.getCenter() ) - atom.getRadius();

        m_cells[cellID]->setEmptyTangentSphere( Sphere(x, y, z, radius) );
        */

        /*
        //  using StringFunctions::strToken()
        rg_INT position = 0;
        string recType  = StringFunctions::subString( currRecLine, position );

        if( recType != QTF_RECORD_VDVTX ) {
            continue;
        }

        rg_INT cellID = atoi( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_REAL x     = atof( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_REAL y     = atof( StringFunctions::strToken( currRecLine, position ).c_str() );
        rg_REAL z     = atof( StringFunctions::strToken( currRecLine, position ).c_str() );

        Sphere  atom  = m_cells[cellID]->getVertex(0)->getBall();

        rg_Point3D center(x, y, z);
        rg_REAL    radius = center.distance( atom.getCenter() ) - atom.getRadius();

        m_cells[cellID]->setEmptyTangentSphere( Sphere(x, y, z, radius) );
        */
    }

}



    
void QTFReader::computeEmptyTangentSpheresOfCells()
{
    m_indexCells.reset4Loop();
    while ( m_indexCells.setNext4Loop() ) {
        QTIndexCell*   currIndexCell = m_indexCells.getpEntity();
        QTTetrahedron* currCell      = currIndexCell->qtCell;

        if ( currCell->isVirtual() || currCell->isIsolatedFace() ) {
            continue;
        }

    
        QTVertex** boundingVertex = currCell->getVertices();
        Sphere     atom[4] = { boundingVertex[0]->getBall(), boundingVertex[1]->getBall(), 
                               boundingVertex[2]->getBall(), boundingVertex[3]->getBall() };

        Sphere emptyTangentSphere[2];
        rg_INT numTS = computeSphereTangentTo4SpheresOutside( atom[0], atom[1], atom[2], atom[3], 
                                                              emptyTangentSphere );

        if ( numTS == 1 ) {
            currCell->setEmptyTangentSphere( emptyTangentSphere[0] );
        }
        else if ( numTS == 2 ) {
            rg_Point3D pointOnTS[4];
            rg_INT i = 0;
            for ( i=0; i<4; i++ ) {
                rg_Point3D vecToAtom = atom[i].getCenter() - emptyTangentSphere[0].getCenter();
                vecToAtom.normalize();

                pointOnTS[i] = emptyTangentSphere[0].getCenter() + (emptyTangentSphere[0].getRadius()*vecToAtom);
            }

            Plane   plane( pointOnTS[1], pointOnTS[2], pointOnTS[3] );
            rg_REAL signedDist = plane.distanceFromPoint( pointOnTS[0] );

            if ( signedDist >= 0.0 ) {
                currCell->setEmptyTangentSphere( emptyTangentSphere[0] );
            }
            else {
                currCell->setEmptyTangentSphere( emptyTangentSphere[1] );
            }
        }
        else {
            //  There is something wrong.
        }
    }
}



rg_BOOL QTFReader::synchronizeAtomsBetweenQTAndMolecule( const Molecule* const& molecule )
{
    rg_BOOL doesSynchronizeQVerticesToAtomsInMolecule = rg_TRUE;

    map<rg_INT, Atom*> atomsInMolecule;

    rg_dList<Atom*> atoms;
    molecule->getPtrAtoms( atoms );


    atoms.reset4Loop();
    while ( atoms.setNext4Loop() ) {
        Atom*  currAtom    = atoms.getEntity();
        rg_INT IDFromInput = currAtom->getSerialFromInputFile();
        atomsInMolecule.insert( make_pair(IDFromInput, currAtom) );
    }


    map <rg_INT, Atom*>::const_iterator it_atom;
    map <rg_INT, Atom*>::const_iterator not_exist = atomsInMolecule.end();


    rg_dList<QTVertex>* vertices = m_QT->getVertices();

    vertices->reset4Loop();
    while(vertices->setNext4Loop()) {
        QTVertex* currQTVertex = vertices->getpEntity();

        if(currQTVertex->isVirtual()) {
            continue;
        }

        Ball*  currBall    = currQTVertex->getBallProperty();
        rg_INT IDFromInput = currBall->getIDFromInput();

        it_atom = atomsInMolecule.find( IDFromInput );
        if ( it_atom != not_exist ) {
            currBall->setProperty( it_atom->second );
        }
        else {
            doesSynchronizeQVerticesToAtomsInMolecule = rg_FALSE;
            break;
        }
    }


    return doesSynchronizeQVerticesToAtomsInMolecule;
}


/*
void QTFReader::extractQuasiTriangulation(list<string>& recordsOfQTFFile, const rg_INT& modelID, list<string*>& recordsOfQTInModel)
{
    rg_BOOL isRecordInSelectedModel = rg_FALSE;

    list<string>::iterator i_recLines = recordsOfQTFFile.begin();
    for( ; i_recLines != recordsOfQTFFile.end(); i_recLines++) {
        string* currRecLine = &(*i_recLines);

        rg_INT position = 0;
        string recType  = StringFunctions::subString( *currRecLine, position );

        if ( !isRecordInSelectedModel ) {
            if( recType == QTF_RECORD_MODEL ) {
                rg_INT  currModelID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );

                if ( modelID == currModelID ) {
                    isRecordInSelectedModel = rg_TRUE;
                }
            }
        }   
        else {
            if( recType == QTF_RECORD_ENDMODEL ) {
                break;
            }

            recordsOfQTInModel.push_back( currRecLine );
        }
    }
}



void QTFReader::readVerticesInQuasiTriangulation(ifstream& fin, list<string>& recordsOfQTFFile)
{
    m_QT->addVertex( QTVertex(0, rg_NULL) );


    list<string>::iterator i_recLines = recordsOfQTFFile.begin();
    while( i_recLines != recordsOfQTFFile.end() ) {
        string* currRecLine = &(*i_recLines);

        rg_INT position = 0;
        string recType  = StringFunctions::subString( *currRecLine, position );

        if( recType != QTF_RECORD_VERTEX ) {
            i_recLines++;
            continue;
        }

        rg_INT  vertexID    = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_INT  firstCellID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() ) - 1;
        rg_REAL x           = atof( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_REAL y           = atof( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_REAL z           = atof( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_REAL radius      = atof( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_INT  IDFromInput = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );

        Ball*     ball   = m_QT->addBall( Ball( vertexID, Sphere(x, y, z, radius), rg_NULL, IDFromInput) );
        QTVertex* vertex = m_QT->addVertex( QTVertex(vertexID, ball) );

        QTIncexVertex indexVertex;
        indexVertex.qtVertex  = vertex;
        indexVertex.firstCell = firstCellID;
        m_indexVertices.add( indexVertex );
        i_recLines++;
    }
}



void QTFReader::readCellsInQuasiTriangulation(ifstream& fin, list<string>& recordsOfQTFFile)
{
    rg_dList< QTIndexCell > qtIndexCellList;

    list<string>::iterator i_recLines = recordsOfQTFFile.begin();

    while( i_recLines != recordsOfQTFFile.end() ) {
        string* currRecLine = &(*i_recLines);

        rg_INT position = 0;
        string recType  = StringFunctions::subString( *currRecLine, position );

        if( recType != QTF_RECORD_CELL ) {
            i_recLines++;
            continue;
        }

        rg_INT  cellID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() ) - 1;

        QTIndexCell indexCell;
        indexCell.vertexIndex[0]   = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        indexCell.vertexIndex[1]   = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        indexCell.vertexIndex[2]   = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        indexCell.vertexIndex[3]   = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        indexCell.neighborIndex[0] = atoi( StringFunctions::subString( *currRecLine, position ).c_str() ) - 1;
        indexCell.neighborIndex[1] = atoi( StringFunctions::subString( *currRecLine, position ).c_str() ) - 1;
        indexCell.neighborIndex[2] = atoi( StringFunctions::subString( *currRecLine, position ).c_str() ) - 1;
        indexCell.neighborIndex[3] = atoi( StringFunctions::subString( *currRecLine, position ).c_str() ) - 1;

        indexCell.qtCell = m_QT->addTetrahedron( QTTetrahedron( cellID ) );
        qtIndexCellList.add( indexCell );

        i_recLines++;
    }



    
    rg_INT numVertices = m_QT->getNumQTVertices() + 1;
    QTVertex** vertex = new QTVertex*[ numVertices ];

    rg_dList<QTVertex>* vertexList = m_QT->getVertices();
    vertexList->reset4Loop();
    while ( vertexList->setNext4Loop() ) {
        QTVertex* currVtx = vertexList->getpEntity();

        vertex[ currVtx->getID() ] = currVtx;
    }


    rg_INT numCells = m_QT->getNumQTTetrahedron();
    QTTetrahedron** cell = new QTTetrahedron*[ numCells ];
    rg_dList<QTTetrahedron>* cellList = m_QT->getTetrahedra();
    cellList->reset4Loop();
    while ( cellList->setNext4Loop() ) {
        QTTetrahedron* currCell = cellList->getpEntity();

        cell[ currCell->getID() ] = currCell;
    }


    m_indexVertices.reset4Loop();
    while ( m_indexVertices.setNext4Loop() ) {
        QTIncexVertex currIndexVtx = m_indexVertices.getEntity();

        currIndexVtx.qtVertex->setFirstTetrahedron( cell[currIndexVtx.firstCell] );
    }


    rg_BOOL doesSetFirstCellOfVirtualVertex = rg_FALSE;

    qtIndexCellList.reset4Loop();
    while ( qtIndexCellList.setNext4Loop() ) {
        QTIndexCell*   currIndexCell = qtIndexCellList.getpEntity();

        QTTetrahedron* currCell = currIndexCell->qtCell;

        if ( currIndexCell->vertexIndex[3] != -1 ) {
            currCell->setVertices(  vertex[ currIndexCell->vertexIndex[0] ], 
                                    vertex[ currIndexCell->vertexIndex[1] ], 
                                    vertex[ currIndexCell->vertexIndex[2] ], 
                                    vertex[ currIndexCell->vertexIndex[3] ] );

            currCell->setNeighbors( cell[  currIndexCell->neighborIndex[0] ], 
                                    cell[  currIndexCell->neighborIndex[1] ], 
                                    cell[  currIndexCell->neighborIndex[2] ], 
                                    cell[  currIndexCell->neighborIndex[3] ] );

            if ( !doesSetFirstCellOfVirtualVertex ) {
                if (    currIndexCell->vertexIndex[0] == 0 
                     || currIndexCell->vertexIndex[1] == 0 
                     || currIndexCell->vertexIndex[2] == 0 
                     || currIndexCell->vertexIndex[3] == 0 ) {
                    vertex[ 0 ]->setFirstTetrahedron( currCell );
                    doesSetFirstCellOfVirtualVertex = rg_TRUE;
                }
            }
//             vertex[ currIndexCell->vertexIndex[0] ]->setFirstTetrahedron( currCell );
//             vertex[ currIndexCell->vertexIndex[1] ]->setFirstTetrahedron( currCell );
//             vertex[ currIndexCell->vertexIndex[2] ]->setFirstTetrahedron( currCell );
//             vertex[ currIndexCell->vertexIndex[3] ]->setFirstTetrahedron( currCell );
        }
        else {
            //  for an isolated triangle.
            currCell->setVertices(  vertex[ currIndexCell->vertexIndex[0] ], 
                                    vertex[ currIndexCell->vertexIndex[1] ], 
                                    vertex[ currIndexCell->vertexIndex[2] ], 
                                    rg_NULL );
            currCell->setNeighbors( rg_NULL, rg_NULL, rg_NULL, rg_NULL );

//             vertex[ currIndexCell->vertexIndex[0] ]->setFirstTetrahedron( currCell );
//             vertex[ currIndexCell->vertexIndex[1] ]->setFirstTetrahedron( currCell );
//             vertex[ currIndexCell->vertexIndex[2] ]->setFirstTetrahedron( currCell );
        }


        if ( currCell->isVirtual() || currCell->isIsolatedFace() ) {
            continue;
        }

        
        QTVertex** boundingVertex = currCell->getVertices();
        Sphere     atom[4] = { boundingVertex[0]->getBall(), boundingVertex[1]->getBall(), 
                               boundingVertex[2]->getBall(), boundingVertex[3]->getBall() };

        Sphere emptyTangentSphere[2];
        rg_INT numTS = computeSphereTangentTo4SpheresOutside( atom[0], atom[1], atom[2], atom[3], 
                                                              emptyTangentSphere );

        if ( numTS == 1 ) {
            currCell->setEmptyTangentSphere( emptyTangentSphere[0] );
        }
        else if ( numTS == 2 ) {
            rg_Point3D pointOnTS[4];
            rg_INT i = 0;
            for ( i=0; i<4; i++ ) {
                rg_Point3D vecToAtom = atom[i].getCenter() - emptyTangentSphere[0].getCenter();
                vecToAtom.normalize();

                pointOnTS[i] = emptyTangentSphere[0].getCenter() + (emptyTangentSphere[0].getRadius()*vecToAtom);
            }

            Plane   plane( pointOnTS[1], pointOnTS[2], pointOnTS[3] );
            rg_REAL signedDist = plane.distanceFromPoint( pointOnTS[0] );

            if ( signedDist >= 0.0 ) {
                currCell->setEmptyTangentSphere( emptyTangentSphere[0] );
            }
            else {
                currCell->setEmptyTangentSphere( emptyTangentSphere[1] );
            }
        }
        else {
            //  There is something wrong.
        }
    }


    delete [] vertex;
    delete [] cell;
}



void QTFReader::readGatesInQuasiTriangulation(ifstream& fin, list<string>& recordsOfQTFFile)
{
    rg_INT numVertices = m_QT->getNumQTVertices() + 1;
    QTVertex** vertex = new QTVertex*[ numVertices ];

    rg_dList<QTVertex>* vertexList = m_QT->getVertices();
    vertexList->reset4Loop();
    while ( vertexList->setNext4Loop() ) {
        QTVertex* currVtx = vertexList->getpEntity();

        vertex[ currVtx->getID() ] = currVtx;
    }


    rg_INT numCells = m_QT->getNumQTTetrahedron();
    QTTetrahedron** cell = new QTTetrahedron*[ numCells ];
    rg_dList<QTTetrahedron>* cellList = m_QT->getTetrahedra();
    cellList->reset4Loop();
    while ( cellList->setNext4Loop() ) {
        QTTetrahedron* currCell = cellList->getpEntity();

        cell[ currCell->getID() ] = currCell;
    }


    
    list<string>::iterator i_recLines = recordsOfQTFFile.begin();

    while( i_recLines != recordsOfQTFFile.end() ) {
        string* currRecLine = &(*i_recLines);

        rg_INT position = 0;
        string recType  = StringFunctions::subString( *currRecLine, position );

        if( recType != QTF_RECORD_GATE ) {
            i_recLines++;
            continue;
        }

        rg_INT length = currRecLine->length();

        rg_INT  gateID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );

        rg_INT vertexID1 = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );
        rg_INT vertexID2 = atoi( StringFunctions::subString( *currRecLine, position ).c_str() );

        QTGate* currGate = m_QT->addGate( QTGate(vertex[vertexID1], vertex[vertexID2]) );

        rg_INT bigCellID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() ) - 1;
        currGate->setBigTetrahedron( cell[bigCellID] );
        
//         vertex[ vertexID1 ]->setFirstTetrahedron( currGate->getBigTetrahedron() );
//         vertex[ vertexID2 ]->setFirstTetrahedron( currGate->getBigTetrahedron() );


        while ( position != -1 ) {
            rg_INT smallWorldID = atoi( StringFunctions::subString( *currRecLine, position ).c_str() ) - 1;
            currGate->addSmallTetrahedron( cell[smallWorldID] );
        }

        i_recLines++;
    }



    delete [] vertex;
    delete [] cell;
}
*/



