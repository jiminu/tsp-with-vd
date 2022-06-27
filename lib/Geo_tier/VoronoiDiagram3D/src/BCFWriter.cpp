#include "BCFWriter.h"
#include "ConstForBetaComplex.h"
#include "StringFunctions.h"
#include "AnalysisOfBetaUniverse.h"
using namespace V::GeometryTier;

#include <time.h>

#include <vector>
#include <algorithm>
#include <map>
using namespace std;

BCFWriter::BCFWriter()
{
}


BCFWriter::~BCFWriter()
{
}
    

void BCFWriter::write(const string& bcfFilename, const string& source, Molecule* molecule, BetaUniverse* QT, const rg_REAL& betaValue, const rg_BOOL& onlyBetaShape)
{
    ofstream fout;
    fout.open( bcfFilename.c_str() );


    rg_INT modelID = 1;

    //  write header of bcf
    modelID = molecule->getModelSerialFromInputFile();
    writeHeaderForPDBMolecule(fout, molecule, betaValue, onlyBetaShape);


    //  write body of bcf
    if ( onlyBetaShape ) {
        writeBetaShape(fout, modelID, QT, betaValue);
    }
    else {
        writeBetaComplex(fout, modelID, QT, betaValue);
    }

    fout.close();
}


void BCFWriter::writeHeaderForPDBMolecule(ofstream& fout, const Molecule* molecule, const rg_REAL& betaValue, const rg_BOOL& onlyBetaShape)
{
    time_t ltime;
    time( &ltime );
    string creatingTime( ctime( &ltime ) );

    fout << BCF_RECORD_SOURCE          << "\t" << "\"" << molecule->getMoleculeFileName()    << "\"" << endl;
    fout << BCF_RECORD_PDBTIMESTAMP    << "\t" << molecule->getTimeStamp()                   << endl;
    fout << BCF_RECORD_PDBFILESIZE     << "\t" << molecule->getFileSize()                    << " BYTE" << endl;
    fout << BCF_RECORD_DEFINE_RADIUS   << "\t" << molecule->getDescriptionOfAtomRadiusType() << endl;

    fout << BCF_RECORD_VDGENENERATOR   << "\t" << "VDRC_VDBALL3D_EDGETRACE_1.0" << endl;
    fout << BCF_RECORD_QTGENENERATOR   << "\t" << "VDRC_QTGEN_1.0" << endl;
    fout << BCF_RECORD_BCFORMATVERSION << "\t" << "1.0" << endl;
    fout << BCF_RECORD_AUTHOR          << "\t" << "VDRC" << "\t" << creatingTime;

    if ( onlyBetaShape ) {
        fout << BCF_RECORD_SHAPE << "\t";
    }
    else {
        fout << BCF_RECORD_COMPLEX << "\t";
    }
    fout << BCF_RECORD_BETAVALUE << " = " << betaValue << endl;
}



void BCFWriter::writeBetaShape(ofstream& fout, const rg_INT& modelID, BetaUniverse* QT, const rg_REAL& betaValue) 
{
    rg_dList<BetaVertex*> verticesInBetaShape;
    rg_dList<BetaEdge*>   edgesInBetaShape;
    rg_dList<BetaFace*>   facesInBetaShape;
	QT->huntVerticesInBetaShape( betaValue, verticesInBetaShape );
	QT->huntEdgesInBetaShape(    betaValue, edgesInBetaShape );
	QT->huntFacesInBetaShape(    betaValue, facesInBetaShape );



    fout << BCF_RECORD_MODEL << "\t" 
         << modelID          << "\t"
         << verticesInBetaShape.getSize()      << "\t"
         << edgesInBetaShape.getSize()         << "\t"
         << facesInBetaShape.getSize()         <<  endl;

    writeStatisticsOfBetaComplexOrShape( fout, QT, betaValue, rg_TRUE);

    map<BetaVertex*, rg_INT> vertexIDs;
    writeVertices(fout, verticesInBetaShape, betaValue, vertexIDs);
    writeEdges(   fout, edgesInBetaShape,    betaValue, vertexIDs);
    writeFaces(   fout, facesInBetaShape,    betaValue, vertexIDs);

    fout << BCF_RECORD_ENDMODEL << endl;
}



void BCFWriter::writeBetaComplex(ofstream& fout, const rg_INT& modelID, BetaUniverse* QT, const rg_REAL& betaValue) 
{
    rg_dList<BetaVertex*> verticesInBetaComplex;
    rg_dList<BetaEdge*>   edgesInBetaComplex;
    rg_dList<BetaFace*>   facesInBetaComplex;
    rg_dList<BetaCell*>   cellsInBetaComplex;
	QT->huntVerticesInBetaComplex( betaValue, verticesInBetaComplex );
	QT->huntEdgesInBetaComplex(    betaValue, edgesInBetaComplex );
	QT->huntFacesInBetaComplex(    betaValue, facesInBetaComplex );
	QT->huntCellsInBetaComplex(    betaValue, cellsInBetaComplex );

    fout << BCF_RECORD_MODEL << "\t" 
         << modelID          << "\t"
         << verticesInBetaComplex.getSize()      << "\t"
         << edgesInBetaComplex.getSize()         << "\t"
         << facesInBetaComplex.getSize()         << "\t"
         << cellsInBetaComplex.getSize()         <<  endl;

    writeStatisticsOfBetaComplexOrShape( fout, QT, betaValue, rg_FALSE);


    map<BetaVertex*, rg_INT> vertexIDs;
    writeVertices(fout, verticesInBetaComplex, betaValue, vertexIDs);
    writeEdges(   fout, edgesInBetaComplex,    betaValue, vertexIDs);
    writeFaces(   fout, facesInBetaComplex,    betaValue, vertexIDs);
    writeCells(   fout, cellsInBetaComplex,    betaValue, vertexIDs);

    fout << BCF_RECORD_ENDMODEL << endl;
}

    

void BCFWriter::writeStatisticsOfBetaComplexOrShape(ofstream& fout, BetaUniverse* QT, const rg_REAL& betaValue, const rg_BOOL& onlyBetaShape)
{
    AnalysisOfBetaUniverse analysisQTBU( QT );

    if ( onlyBetaShape ) {
        rg_INT numVtxsInBetaComplex[4]  = {0, 0, 0, 0};
        rg_INT numEdgesInBetaComplex[4] = {0, 0, 0, 0};
        rg_INT numFacesInBetaComplex[4] = {0, 0, 0, 0};
        analysisQTBU.countVerticesInBetaComplex( betaValue, numVtxsInBetaComplex[0],  numVtxsInBetaComplex[1],  numVtxsInBetaComplex[2],  numVtxsInBetaComplex[3]);
        analysisQTBU.countEdgesInBetaComplex(    betaValue, numEdgesInBetaComplex[0], numEdgesInBetaComplex[1], numEdgesInBetaComplex[2], numEdgesInBetaComplex[3]);
        analysisQTBU.countFacesInBetaComplex(    betaValue, numFacesInBetaComplex[0], numFacesInBetaComplex[1], numFacesInBetaComplex[2], numFacesInBetaComplex[3]);


        rg_INT numVtxInBetaShape   = numVtxsInBetaComplex[1]  + numVtxsInBetaComplex[2];
        rg_INT numEdgesInBetaShape = numEdgesInBetaComplex[1] + numEdgesInBetaComplex[2];
        rg_INT numFacesInBetaShape = numFacesInBetaComplex[1] + numFacesInBetaComplex[2];

        fout << BCF_RECORD_REMARK << "\t" << "# vertices = " << numVtxInBetaShape   << "\t" << "(S: " << numVtxsInBetaComplex[1]  << ", R: " << numVtxsInBetaComplex[2]  << ")" << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# edges    = " << numEdgesInBetaShape << "\t" << "(S: " << numEdgesInBetaComplex[1] << ", R: " << numEdgesInBetaComplex[2] << ")" << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# faces    = " << numFacesInBetaShape << "\t" << "(S: " << numFacesInBetaComplex[1] << ", R: " << numFacesInBetaComplex[2] << ")" << endl;
    }
    else {
        rg_INT numVtxsInBetaComplex[4]  = {0, 0, 0, 0};
        rg_INT numEdgesInBetaComplex[4] = {0, 0, 0, 0};
        rg_INT numFacesInBetaComplex[4] = {0, 0, 0, 0};
        rg_INT numCellsInBetaComplex[2] = {0, 0};
        analysisQTBU.countVerticesInBetaComplex( betaValue, numVtxsInBetaComplex[0],  numVtxsInBetaComplex[1],  numVtxsInBetaComplex[2],  numVtxsInBetaComplex[3]);
        analysisQTBU.countEdgesInBetaComplex(    betaValue, numEdgesInBetaComplex[0], numEdgesInBetaComplex[1], numEdgesInBetaComplex[2], numEdgesInBetaComplex[3]);
        analysisQTBU.countFacesInBetaComplex(    betaValue, numFacesInBetaComplex[0], numFacesInBetaComplex[1], numFacesInBetaComplex[2], numFacesInBetaComplex[3]);
        analysisQTBU.countCellsInBetaComplex(    betaValue, numCellsInBetaComplex[0], numCellsInBetaComplex[1]);


        //////////////////////////////////////////////////////////////
        //  for # simplexes
        rg_INT numAllVtxInBetaComplex   = numVtxsInBetaComplex[1]  + numVtxsInBetaComplex[2]  + numVtxsInBetaComplex[3];
        rg_INT numAllEdgesInBetaComplex = numEdgesInBetaComplex[1] + numEdgesInBetaComplex[2] + numEdgesInBetaComplex[3];
        rg_INT numAllFacesInBetaComplex = numFacesInBetaComplex[1] + numFacesInBetaComplex[2] + numFacesInBetaComplex[3];
        fout << BCF_RECORD_REMARK << "\t" << "# vertices = " << numAllVtxInBetaComplex   << "\t" << "(S: " << numVtxsInBetaComplex[1]  << ", R: " << numVtxsInBetaComplex[2]  << ", I: " << numVtxsInBetaComplex[3] << ")" << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# edges    = " << numAllEdgesInBetaComplex << "\t" << "(S: " << numEdgesInBetaComplex[1] << ", R: " << numEdgesInBetaComplex[2] << ", I: " << numEdgesInBetaComplex[3] << ")" << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# faces    = " << numAllFacesInBetaComplex << "\t" << "(S: " << numFacesInBetaComplex[1] << ", R: " << numFacesInBetaComplex[2] << ", I: " << numFacesInBetaComplex[3] << ")" << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# cells    = " << numCellsInBetaComplex[1] << "\t" << "(I: " << numCellsInBetaComplex[1] << ")" << endl;

        rg_INT numCellsWithMinusVolumeInBetaComplex[2] = {0, 0};
        analysisQTBU.countCellsWithMinusVolumeInBetaComplex( betaValue, numCellsWithMinusVolumeInBetaComplex[0], numCellsWithMinusVolumeInBetaComplex[1]);
        fout << BCF_RECORD_REMARK << "\t" << "(# cells with a minus volume = " << numCellsWithMinusVolumeInBetaComplex[1] << ")" << endl;


        //////////////////////////////////////////////////////////////
        //  for anomaly
        rg_INT num2AdjacencyAnomaly[3]  = {0, 0, 0};
        rg_INT num3AdjacencyAnomaly[3]  = {0, 0, 0};
        rg_INT num4AdjacencyAnomaly[3]  = {0, 0, 0};
        analysisQTBU.countCellStateOf2AdjacencyAnomaliesInBetaComplex(betaValue, num2AdjacencyAnomaly[0], num2AdjacencyAnomaly[1], num2AdjacencyAnomaly[2]);
        analysisQTBU.countCellStateOf3AdjacencyAnomaliesInBetaComplex(betaValue, num3AdjacencyAnomaly[0], num3AdjacencyAnomaly[1], num3AdjacencyAnomaly[2]);
        analysisQTBU.countCellStateOf4AdjacencyAnomaliesInBetaComplex(betaValue, num4AdjacencyAnomaly[0], num4AdjacencyAnomaly[1], num4AdjacencyAnomaly[2]);
        fout << BCF_RECORD_REMARK << "\t" << "# 2-adjacency anomaly = " << num2AdjacencyAnomaly[2] << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# 3-adjacency anomaly = " << num3AdjacencyAnomaly[2] << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# 4-adjacency anomaly = " << num4AdjacencyAnomaly[2] << endl;


        rg_INT numDanglingAnomaly[3]  = {0, 0, 0};
        analysisQTBU.countDanglingAnomaliesInBetaComplex( betaValue, numDanglingAnomaly[0], numDanglingAnomaly[1], numDanglingAnomaly[2]);

        fout << BCF_RECORD_REMARK << "\t" << "# dangling-face anomaly    = " << numDanglingAnomaly[0] << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# dangling-cell anomaly    = " << numDanglingAnomaly[1] << endl;
        fout << BCF_RECORD_REMARK << "\t" << "# dangling-cluster anomaly = " << numDanglingAnomaly[2] << endl;

    }
}



void BCFWriter::writeVertices(ofstream& fout, const rg_dList<BetaVertex*>& vertices, const rg_REAL& betaValue, map<BetaVertex*, rg_INT>& vertexIDs ) 
{
    vector<BetaVertex*> vertexList;
    vertices.reset4Loop();
    while ( vertices.setNext4Loop() ) {
        BetaVertex* currVtx = vertices.getEntity();
        vertexList.push_back( currVtx );
    }
    sort( vertexList.begin(), vertexList.end(), BCVertexInputIDLess );



    rg_INT vertexID = 1;
    vector<BetaVertex*>::iterator i_vtx;
    for ( i_vtx=vertexList.begin(); i_vtx!=vertexList.end(); i_vtx++, vertexID++ ) {
        BetaVertex* currVtx = *i_vtx;

        rg_Point3D center = currVtx->getBall().getCenter();
        rg_REAL    radius = currVtx->getBall().getRadius();

        fout << BCF_RECORD_VERTEX            << "\t";
        fout << vertexID  << "\t"; 
        fout << center.getX()              << "\t";
        fout << center.getY()              << "\t";
        fout << center.getZ()              << "\t";         
        fout << radius                     << "\t"; 

        rg_INT state = currVtx->getBoundingState(betaValue);
        switch ( state ) {
        case SINGULAR_SIMPLEX:
            fout << BCF_SINGULAR_STATE << "\t";
            break;
        case REGULAR_SIMPLEX:
            fout << BCF_REGULAR_STATE << "\t";
            break;
        case INTERIOR_SIMPLEX:
            fout << BCF_INTERIOR_STATE << "\t";
            break;
        }
        fout << currVtx->getBallProperty()->getIDFromInput() << endl;

        vertexIDs.insert( make_pair(currVtx, vertexID) );
    }
}



void BCFWriter::writeEdges(ofstream& fout, const rg_dList<BetaEdge*>& edges, const rg_REAL& betaValue, const map<BetaVertex*, rg_INT>& vertexIDs ) 
{
    vector<BetaEdgeByVertexID> edgeArray;
    edges.reset4Loop();
    while ( edges.setNext4Loop() ) {
        BetaEdge* currEdge = edges.getEntity();

        rg_INT startVtxID = vertexIDs.find( currEdge->getStartVertex() )->second;
        rg_INT endVtxID   = vertexIDs.find( currEdge->getEndVertex() )->second;
        edgeArray.push_back( BetaEdgeByVertexID(currEdge, startVtxID, endVtxID ) );
    }
    sort( edgeArray.begin(), edgeArray.end() );


    rg_INT edgeID = 1;
    vector<BetaEdgeByVertexID>::iterator i_edge;
    for ( i_edge=edgeArray.begin(); i_edge!=edgeArray.end(); i_edge++, edgeID++ ) {
        BetaEdgeByVertexID* currEdge = &(*i_edge);

        fout << BCF_RECORD_EDGE            << "\t";
        fout << edgeID  << "\t"; 
        fout << currEdge->getVertexID(0)   << "\t";
        fout << currEdge->getVertexID(1)   << "\t";

        rg_INT state = currEdge->getEdge()->getBoundingState(betaValue);
        switch ( state ) {
        case SINGULAR_SIMPLEX:
            fout << BCF_SINGULAR_STATE;
            break;
        case REGULAR_SIMPLEX:
            fout << BCF_REGULAR_STATE;
            break;
        case INTERIOR_SIMPLEX:
            fout << BCF_INTERIOR_STATE;
            break;
        }
        fout << endl;
    }
}



void BCFWriter::writeFaces(ofstream& fout, const rg_dList<BetaFace*>& faces, const rg_REAL& betaValue, const map<BetaVertex*, rg_INT>& vertexIDs ) 
{
    vector<BetaFaceByVertexID> faceArray;
    faces.reset4Loop();
    while ( faces.setNext4Loop() ) {
        BetaFace* currFace = faces.getEntity();

        BetaVertex* vertex[3] = {rg_NULL, rg_NULL, rg_NULL};
        currFace->searchVerticesInIntraWorld( vertex );

        rg_INT vtxID[3];
        for ( rg_INT i=0; i<3; ++i ) {
            vtxID[i] = vertexIDs.find( vertex[i] )->second;
        }

        if ( currFace->getBoundingState(betaValue) == REGULAR_SIMPLEX ) {
            if ( currFace->getLeftCell()->getBoundingState(betaValue) == INTERIOR_SIMPLEX ) {
                faceArray.push_back( BetaFaceByVertexID(currFace, vtxID[0], vtxID[1], vtxID[2] ) );
            }
            else {
                faceArray.push_back( BetaFaceByVertexID(currFace, vtxID[2], vtxID[1], vtxID[0] ) );
            }
        }
        else {
            faceArray.push_back( BetaFaceByVertexID(currFace, vtxID[0], vtxID[1], vtxID[2] ) );
        }
    }
    sort( faceArray.begin(), faceArray.end() );


    rg_INT faceID = 1;
    vector<BetaFaceByVertexID>::iterator i_face;
    for ( i_face=faceArray.begin(); i_face!=faceArray.end(); i_face++, faceID++ ) {
        BetaFaceByVertexID* currFace = &(*i_face);

        fout << BCF_RECORD_FACE            << "\t";
        fout << faceID  << "\t"; 
        fout << currFace->getVertexID(0)   << "\t";
        fout << currFace->getVertexID(1)   << "\t";
        fout << currFace->getVertexID(2)   << "\t";

        rg_INT state = currFace->getFace()->getBoundingState(betaValue);
        switch ( state ) {
        case SINGULAR_SIMPLEX:
            fout << BCF_SINGULAR_STATE;
            break;
        case REGULAR_SIMPLEX:
            fout << BCF_REGULAR_STATE;
            break;
        case INTERIOR_SIMPLEX:
            fout << BCF_INTERIOR_STATE;
            break;
        }
        fout << endl;
    }
}



void BCFWriter::writeCells(ofstream& fout, const rg_dList<BetaCell*>& cells, const rg_REAL& betaValue, const map<BetaVertex*, rg_INT>& vertexIDs ) 
{
    vector<BetaCellByVertexID> cellArray;
    cells.reset4Loop();
    while ( cells.setNext4Loop() ) {
        BetaCell* currCell = cells.getEntity();

        rg_INT vtxID[4];
        for ( rg_INT i=0; i<4; ++i ) {
            vtxID[i] = vertexIDs.find( currCell->getVertex(i) )->second;
        }

        bool volumeSign = true;
        if ( rg_NEG( currCell->computeSignedVolume() ) ) {
            volumeSign = false;
        }

        cellArray.push_back( BetaCellByVertexID(currCell, vtxID[0], vtxID[1], vtxID[2], vtxID[3], volumeSign ) );
    }
    sort( cellArray.begin(), cellArray.end() );


    rg_INT cellID = 1;
    vector<BetaCellByVertexID>::iterator i_cell;
    for ( i_cell=cellArray.begin(); i_cell!=cellArray.end(); i_cell++, cellID++ ) {
        BetaCellByVertexID* currCell = &(*i_cell);

        fout << BCF_RECORD_CELL            << "\t";
        fout << cellID  << "\t"; 
        fout << currCell->getVertexID(0)   << "\t";
        fout << currCell->getVertexID(1)   << "\t";
        fout << currCell->getVertexID(2)   << "\t";
        fout << currCell->getVertexID(3)   << "\t";

        fout << BCF_INTERIOR_STATE;

        if ( !currCell->hasPositiveVolume() ) {
            fout << "\t*";
        }
        fout << endl;
    }
}






/////////////////////////////////////////////////////////////////////////////////

BCFWriter::BetaEdgeByVertexID::BetaEdgeByVertexID()
{
    m_edge      = rg_NULL;
    m_vtxIDs[0] = -1;;
    m_vtxIDs[1] = -1;;
}



BCFWriter::BetaEdgeByVertexID::BetaEdgeByVertexID(BetaEdge* edge, const rg_INT& vtxID1, const rg_INT& vtxID2)
{
    m_edge      = edge;     
    m_vtxIDs[0] = vtxID1;
    m_vtxIDs[1] = vtxID2;
    if ( vtxID1 > vtxID2 ) {
        m_vtxIDs[0] = vtxID2;
        m_vtxIDs[1] = vtxID1;
    }
}



BCFWriter::BetaEdgeByVertexID::BetaEdgeByVertexID(const BetaEdgeByVertexID& betaEdge)
{
    m_edge      = betaEdge.m_edge;     
    m_vtxIDs[0] = betaEdge.m_vtxIDs[0];
    m_vtxIDs[1] = betaEdge.m_vtxIDs[1];
}



BCFWriter::BetaEdgeByVertexID::~BetaEdgeByVertexID()
{
}




void      BCFWriter::BetaEdgeByVertexID::setBetaEdgeByVertexID(BetaEdge* edge, const rg_INT& vtxID1, const rg_INT& vtxID2)
{
    m_edge      = edge;     
    m_vtxIDs[0] = vtxID1;
    m_vtxIDs[1] = vtxID2;
    if ( vtxID1 > vtxID2 ) {
        m_vtxIDs[0] = vtxID2;
        m_vtxIDs[1] = vtxID1;
    }
}



BCFWriter::BetaEdgeByVertexID& BCFWriter::BetaEdgeByVertexID::operator =(const BetaEdgeByVertexID& betaEdge)
{
    if ( this == &betaEdge ) {
        return *this;
    }

    m_edge      = betaEdge.m_edge;     
    m_vtxIDs[0] = betaEdge.m_vtxIDs[0];
    m_vtxIDs[1] = betaEdge.m_vtxIDs[1];

    return *this;
}



bool BCFWriter::BetaEdgeByVertexID::operator==(const BetaEdgeByVertexID& betaEdge) const
{
    if ( m_edge == betaEdge.m_edge && m_vtxIDs[0] == betaEdge.m_vtxIDs[0] && m_vtxIDs[1] == betaEdge.m_vtxIDs[1] ) {
        return true;
    }
    else {
        return false;
    }
}



bool BCFWriter::BetaEdgeByVertexID::operator <(const BetaEdgeByVertexID& betaEdge) const
{
    if ( m_vtxIDs[0] < betaEdge.m_vtxIDs[0] ) {
        return true;
    }
    else if ( m_vtxIDs[0] == betaEdge.m_vtxIDs[0] ) {
        if ( m_vtxIDs[1] < betaEdge.m_vtxIDs[1] ) {
            return true;
        }
    }

    return false;
}



BCFWriter::BetaFaceByVertexID::BetaFaceByVertexID()
{
    m_face      = rg_NULL;
    m_vtxIDs[0] = -1;
    m_vtxIDs[1] = -1;
    m_vtxIDs[2] = -1;
}



BCFWriter::BetaFaceByVertexID::BetaFaceByVertexID(BetaFace* face, const rg_INT& vtxID1, const rg_INT& vtxID2, const rg_INT& vtxID3)
{
    m_face      = face;
    if ( vtxID1 <= vtxID2 && vtxID1 <= vtxID3 ) {
        m_vtxIDs[0] = vtxID1;
        m_vtxIDs[1] = vtxID2;
        m_vtxIDs[2] = vtxID3;
    }
    else if ( vtxID2 <= vtxID1 && vtxID2 <= vtxID3 ) {
        m_vtxIDs[0] = vtxID2;
        m_vtxIDs[1] = vtxID3;
        m_vtxIDs[2] = vtxID1;
    }
    else { // if ( vtxID3 <= vtxID1 && vtxID3 <= vtxID2 ) {
        m_vtxIDs[0] = vtxID3;
        m_vtxIDs[1] = vtxID1;
        m_vtxIDs[2] = vtxID2;
    }
}



BCFWriter::BetaFaceByVertexID::BetaFaceByVertexID(const BetaFaceByVertexID& betaFace)
{
    m_face      = betaFace.m_face;
    m_vtxIDs[0] = betaFace.m_vtxIDs[0];
    m_vtxIDs[1] = betaFace.m_vtxIDs[1];
    m_vtxIDs[2] = betaFace.m_vtxIDs[2];
}



BCFWriter::BetaFaceByVertexID::~BetaFaceByVertexID()
{
}



void      BCFWriter::BetaFaceByVertexID::setBetaFaceByVertexID(BetaFace* face, const rg_INT& vtxID1, const rg_INT& vtxID2, const rg_INT& vtxID3)
{
    m_face      = face;
    if ( vtxID1 <= vtxID2 && vtxID1 <= vtxID3 ) {
        m_vtxIDs[0] = vtxID1;
        m_vtxIDs[1] = vtxID2;
        m_vtxIDs[2] = vtxID3;
    }
    else if ( vtxID2 <= vtxID1 && vtxID2 <= vtxID3 ) {
        m_vtxIDs[0] = vtxID2;
        m_vtxIDs[1] = vtxID3;
        m_vtxIDs[2] = vtxID1;
    }
    else { // if ( vtxID3 <= vtxID1 && vtxID3 <= vtxID2 ) {
        m_vtxIDs[0] = vtxID3;
        m_vtxIDs[1] = vtxID1;
        m_vtxIDs[2] = vtxID2;
    }
}



BCFWriter::BetaFaceByVertexID& BCFWriter::BetaFaceByVertexID::operator =(const BetaFaceByVertexID& betaFace)
{
    if ( this == &betaFace ) {
        return *this;
    }

    m_face      = betaFace.m_face;
    m_vtxIDs[0] = betaFace.m_vtxIDs[0];
    m_vtxIDs[1] = betaFace.m_vtxIDs[1];
    m_vtxIDs[2] = betaFace.m_vtxIDs[2];

    return *this;
}



bool BCFWriter::BetaFaceByVertexID::operator==(const BetaFaceByVertexID& betaFace) const
{
    if ( m_face == betaFace.m_face && m_vtxIDs[0] == betaFace.m_vtxIDs[0] && m_vtxIDs[1] == betaFace.m_vtxIDs[1] && m_vtxIDs[2] == betaFace.m_vtxIDs[2] ) {
        return true;
    }
    else {
        return false;
    }
}



bool BCFWriter::BetaFaceByVertexID::operator <(const BetaFaceByVertexID& betaFace) const
{
    if ( m_vtxIDs[0] < betaFace.m_vtxIDs[0] ) {
        return true;
    }
    else if ( m_vtxIDs[0] == betaFace.m_vtxIDs[0] ) {
        if ( m_vtxIDs[1] < betaFace.m_vtxIDs[1] ) {
            return true;
        }
        else if ( m_vtxIDs[1] == betaFace.m_vtxIDs[1] ) {
            if ( m_vtxIDs[2] < betaFace.m_vtxIDs[2] ) {
                return true;
            }
        }
    }

    return false;
}




BCFWriter::BetaCellByVertexID::BetaCellByVertexID()
{
    m_cell      = rg_NULL;
    m_vtxIDs[0] = -1;
    m_vtxIDs[1] = -1;
    m_vtxIDs[2] = -1;
    m_vtxIDs[3] = -1;

    m_positiveVolume = true;
}



BCFWriter::BetaCellByVertexID::BetaCellByVertexID(BetaCell* cell, const rg_INT& vtxID1, const rg_INT& vtxID2, const rg_INT& vtxID3, const rg_INT& vtxID4, const bool& positiveVolume)
{
    m_cell      = cell;
    m_positiveVolume = positiveVolume;

    if ( vtxID1 <= vtxID2 && vtxID1 <= vtxID3 && vtxID1 <= vtxID4) {
        m_vtxIDs[0] = vtxID1;

        if ( vtxID2 <= vtxID3 && vtxID2 <= vtxID4) {
            m_vtxIDs[1] = vtxID2;
            m_vtxIDs[2] = vtxID3;
            m_vtxIDs[3] = vtxID4;
        }
        else if ( vtxID3 <= vtxID2 && vtxID3 <= vtxID4) {
            m_vtxIDs[1] = vtxID3;
            m_vtxIDs[2] = vtxID4;
            m_vtxIDs[3] = vtxID2;
        }
        else {
            m_vtxIDs[1] = vtxID4;
            m_vtxIDs[2] = vtxID2;
            m_vtxIDs[3] = vtxID3;
        }
    }
    else if ( vtxID2 <= vtxID1 && vtxID2 <= vtxID3 && vtxID2 <= vtxID4) {
        m_vtxIDs[0] = vtxID2;
            
        if ( vtxID1 <= vtxID3 && vtxID1 <= vtxID4) {
            m_vtxIDs[1] = vtxID1;
            m_vtxIDs[2] = vtxID4;
            m_vtxIDs[3] = vtxID3;
        }
        else if ( vtxID4 <= vtxID1 && vtxID4 <= vtxID3) {
            m_vtxIDs[1] = vtxID4;
            m_vtxIDs[2] = vtxID3;
            m_vtxIDs[3] = vtxID1;
        }
        else {
            m_vtxIDs[1] = vtxID3;
            m_vtxIDs[2] = vtxID1;
            m_vtxIDs[3] = vtxID4;
        }
    }
    else if ( vtxID3 <= vtxID1 && vtxID3 <= vtxID2 && vtxID3 <= vtxID4) {
        m_vtxIDs[0] = vtxID3;

        if ( vtxID1 <= vtxID2 && vtxID1 <= vtxID4) {
            m_vtxIDs[1] = vtxID1;
            m_vtxIDs[2] = vtxID2;
            m_vtxIDs[3] = vtxID4;
        }
        else if ( vtxID2 <= vtxID1 && vtxID2 <= vtxID4) {
            m_vtxIDs[1] = vtxID2;
            m_vtxIDs[2] = vtxID4;
            m_vtxIDs[3] = vtxID1;
        }
        else {
            m_vtxIDs[1] = vtxID4;
            m_vtxIDs[2] = vtxID1;
            m_vtxIDs[3] = vtxID2;
        }
    }
    else { //if ( vtxID4 <= vtxID1 && vtxID4 <= vtxID2 && vtxID4 <= vtxID3) {
        m_vtxIDs[0] = vtxID4;

        if ( vtxID1 <= vtxID2 && vtxID1 <= vtxID3) {
            m_vtxIDs[1] = vtxID1;
            m_vtxIDs[2] = vtxID3;
            m_vtxIDs[3] = vtxID2;
        }
        else if ( vtxID3 <= vtxID1 && vtxID3 <= vtxID2) {
            m_vtxIDs[1] = vtxID3;
            m_vtxIDs[2] = vtxID2;
            m_vtxIDs[3] = vtxID1;
        }
        else {
            m_vtxIDs[1] = vtxID2;
            m_vtxIDs[2] = vtxID1;
            m_vtxIDs[3] = vtxID3;
        }
    }
}



BCFWriter::BetaCellByVertexID::BetaCellByVertexID(const BetaCellByVertexID& betaCell)
{
    m_cell      = betaCell.m_cell;
    m_vtxIDs[0] = betaCell.m_vtxIDs[0];
    m_vtxIDs[1] = betaCell.m_vtxIDs[1];
    m_vtxIDs[2] = betaCell.m_vtxIDs[2];
    m_vtxIDs[3] = betaCell.m_vtxIDs[3];
    m_positiveVolume = betaCell.m_positiveVolume;
}



BCFWriter::BetaCellByVertexID::~BetaCellByVertexID()
{
}



void      BCFWriter::BetaCellByVertexID::setBetaCellByVertexID(BetaCell* cell, const rg_INT& vtxID1, const rg_INT& vtxID2, const rg_INT& vtxID3, const rg_INT& vtxID4, const bool& positiveVolume)
{
    m_cell      = cell;
    m_positiveVolume = positiveVolume;

    if ( vtxID1 <= vtxID2 && vtxID1 <= vtxID3 && vtxID1 <= vtxID4) {
        m_vtxIDs[0] = vtxID1;

        if ( vtxID2 <= vtxID3 && vtxID2 <= vtxID4) {
            m_vtxIDs[1] = vtxID2;
            m_vtxIDs[2] = vtxID3;
            m_vtxIDs[3] = vtxID4;
        }
        else if ( vtxID3 <= vtxID2 && vtxID3 <= vtxID4) {
            m_vtxIDs[1] = vtxID3;
            m_vtxIDs[2] = vtxID4;
            m_vtxIDs[3] = vtxID2;
        }
        else {
            m_vtxIDs[1] = vtxID4;
            m_vtxIDs[2] = vtxID2;
            m_vtxIDs[3] = vtxID3;
        }
    }
    else if ( vtxID2 <= vtxID1 && vtxID2 <= vtxID3 && vtxID2 <= vtxID4) {
        m_vtxIDs[0] = vtxID2;
            
        if ( vtxID1 <= vtxID3 && vtxID1 <= vtxID4) {
            m_vtxIDs[1] = vtxID1;
            m_vtxIDs[2] = vtxID4;
            m_vtxIDs[3] = vtxID3;
        }
        else if ( vtxID4 <= vtxID1 && vtxID4 <= vtxID3) {
            m_vtxIDs[1] = vtxID4;
            m_vtxIDs[2] = vtxID3;
            m_vtxIDs[3] = vtxID1;
        }
        else {
            m_vtxIDs[1] = vtxID3;
            m_vtxIDs[2] = vtxID1;
            m_vtxIDs[3] = vtxID4;
        }
    }
    else if ( vtxID3 <= vtxID1 && vtxID3 <= vtxID2 && vtxID3 <= vtxID4) {
        m_vtxIDs[0] = vtxID3;

        if ( vtxID1 <= vtxID2 && vtxID1 <= vtxID4) {
            m_vtxIDs[1] = vtxID1;
            m_vtxIDs[2] = vtxID2;
            m_vtxIDs[3] = vtxID4;
        }
        else if ( vtxID2 <= vtxID1 && vtxID2 <= vtxID4) {
            m_vtxIDs[1] = vtxID2;
            m_vtxIDs[2] = vtxID4;
            m_vtxIDs[3] = vtxID1;
        }
        else {
            m_vtxIDs[1] = vtxID4;
            m_vtxIDs[2] = vtxID1;
            m_vtxIDs[3] = vtxID2;
        }
    }
    else { //if ( vtxID4 <= vtxID1 && vtxID4 <= vtxID2 && vtxID4 <= vtxID3) {
        m_vtxIDs[0] = vtxID4;

        if ( vtxID1 <= vtxID2 && vtxID1 <= vtxID3) {
            m_vtxIDs[1] = vtxID1;
            m_vtxIDs[2] = vtxID3;
            m_vtxIDs[3] = vtxID2;
        }
        else if ( vtxID3 <= vtxID1 && vtxID3 <= vtxID2) {
            m_vtxIDs[1] = vtxID3;
            m_vtxIDs[2] = vtxID2;
            m_vtxIDs[3] = vtxID1;
        }
        else {
            m_vtxIDs[1] = vtxID2;
            m_vtxIDs[2] = vtxID1;
            m_vtxIDs[3] = vtxID3;
        }
    }
}



BCFWriter::BetaCellByVertexID& BCFWriter::BetaCellByVertexID::operator =(const BetaCellByVertexID& betaCell)
{
    if ( this == &betaCell ) {
        return *this;
    }

    m_cell      = betaCell.m_cell;
    m_vtxIDs[0] = betaCell.m_vtxIDs[0];
    m_vtxIDs[1] = betaCell.m_vtxIDs[1];
    m_vtxIDs[2] = betaCell.m_vtxIDs[2];
    m_vtxIDs[3] = betaCell.m_vtxIDs[3];
    m_positiveVolume = betaCell.m_positiveVolume;

    return *this;
}



bool BCFWriter::BetaCellByVertexID::operator==(const BetaCellByVertexID& betaCell) const
{
    if ( m_cell == betaCell.m_cell && m_vtxIDs[0] == betaCell.m_vtxIDs[0] && m_vtxIDs[1] == betaCell.m_vtxIDs[1] && m_vtxIDs[2] == betaCell.m_vtxIDs[2] && m_vtxIDs[3] == betaCell.m_vtxIDs[3]) {
        return true;
    }
    else {
        return false;
    }
}



bool BCFWriter::BetaCellByVertexID::operator <(const BetaCellByVertexID& betaCell) const
{
    if ( m_vtxIDs[0] < betaCell.m_vtxIDs[0] ) {
        return true;
    }
    else if ( m_vtxIDs[0] == betaCell.m_vtxIDs[0] ) {
        if ( m_vtxIDs[1] < betaCell.m_vtxIDs[1] ) {
            return true;
        }
        else if ( m_vtxIDs[1] == betaCell.m_vtxIDs[1] ) {
            if ( m_vtxIDs[2] < betaCell.m_vtxIDs[2] ) {
                return true;
            }
            else if ( m_vtxIDs[2] == betaCell.m_vtxIDs[2] ) {
                if ( m_vtxIDs[3] < betaCell.m_vtxIDs[3] ) {
                    return true;
                }
            }
        }
    }

    return false;
}



