#include <time.h>
#include "rg_QTFWriter.h"
#include "ConstForQuasiTriangulation.h"
#include "StringFunctions.h"
using namespace V::GeometryTier;


QTFWriter::QTFWriter()
: m_QT(rg_NULL)
{
}



QTFWriter::QTFWriter(QuasiTriangulation* QT)
: m_QT(QT)
{
}



QTFWriter::QTFWriter(const QTFWriter& writer)
: m_QT(writer.m_QT)
{
}



QTFWriter::~QTFWriter()
{
    m_QT = rg_NULL;
}




QuasiTriangulation* QTFWriter::getQuasiTriangulation() const
{   
    return m_QT;
}



void QTFWriter::setQuasiTriangulation(QuasiTriangulation* QT)
{
    m_QT = QT;
}




QTFWriter& QTFWriter::operator=(const QTFWriter& writer)
{
    if ( this == &writer) {
        return *this;
    }

    m_QT = writer.m_QT;

    return *this;
}



void QTFWriter::write( const string& qtfFilename, const string& source, QuasiTriangulation* QT, Molecule* molecule )
{
    QT->resetIDs();


    ofstream fout;

    if ( m_QTFileName.empty() ) {
        m_QTFileName = qtfFilename;
        
        fout.open( m_QTFileName.c_str() );
        
        rg_INT modelID = 0;
        if ( molecule == rg_NULL ) {
            modelID = 1;
            writeHeaderForMolecule(fout, source, "" );
        }
        else {
            modelID = molecule->getModelSerialFromInputFile();

            string moleculeFilename      = molecule->getMoleculeFileName();
            string moleculeFileExtension = StringFunctions::strToLower( StringFunctions::getExtensionFromFileName( moleculeFilename ) );   

            if ( moleculeFileExtension == "pdb" || moleculeFileExtension == "ent") {
                writeHeaderForPDBMolecule(fout, molecule);
            }
            else {
                writeHeaderForMolecule(fout, source, molecule->getDescriptionOfAtomRadiusType() );
            }
        }

        writeQuasiTriangulation(fout, modelID, QT);

        fout.close();
    }
    else {
        if ( m_QTFileName != qtfFilename ) {
            return;
        }

        ofstream fout;
        fout.open( m_QTFileName.c_str(), ios_base::out|ios_base::app );
         
        rg_INT modelID = 0;
        if ( molecule == rg_NULL ) {
            modelID = 1;
        }
        else {
            modelID = molecule->getModelSerialFromInputFile();
        }

        writeQuasiTriangulation(fout, modelID, QT);

        fout.close();
    }
}



void QTFWriter::write( const string& qtfFilename, 
                       const string& source,
                       const rg_dList< pair<QuasiTriangulation*, Molecule*> >& models )
{


    pair<QuasiTriangulation*, Molecule*>* firstModel = rg_NULL;

    models.reset4Loop();
    while ( models.setNext4Loop() ) {
        pair<QuasiTriangulation*, Molecule*>* currModel = models.getpEntity();
        if ( currModel->first == rg_NULL || currModel->second == rg_NULL ) {
            continue;
        }

        if ( firstModel == rg_NULL ) {
            firstModel = currModel;
        }
    }

    if ( firstModel == rg_NULL ) {
        return;
    }



    ofstream fout;
    fout.open( qtfFilename.c_str() );

    writeHeaderForMolecule(fout, source, firstModel->second->getDescriptionOfAtomRadiusType());

    models.reset4Loop();
    while ( models.setNext4Loop() ) {
        pair<QuasiTriangulation*, Molecule*> currModel = models.getEntity();

        if ( currModel.first == rg_NULL || currModel.second == rg_NULL ) {
            continue;
        }


        currModel.first->resetIDs();


        rg_INT modelID = currModel.second->getModelSerialFromInputFile();

        writeQuasiTriangulation(fout, modelID, currModel.first);
    }

    fout.close();
}


/*
void QTFWriter::write(const string& qtfFilename, const string& source)
{
    m_QT->resetIDs();

    ofstream fout;
    fout.open( qtfFilename.c_str() );

    writeSource(fout, source);
    writeVerticesInQuasiTriangulation(fout);
    writeCellsInQuasiTriangulation(fout);
    writeGatesInQuasiTriangulation(fout);

    fout.close();
}
*/


rg_INT QTFWriter::getVertexID(QuasiTriangulation* QT, QTVertex* vertex) 
{
    rg_INT vertexID = -1;
    
    if ( vertex != rg_NULL ) {
        vertexID = QT->getAbsoluteIDOfVertex( vertex );
    }

    return vertexID;
}



rg_INT QTFWriter::getCellID(QTTetrahedron* cell) 
{
    rg_INT cellID = cell->getID() + 1;

    return cellID;
}


/*
void QTFWriter::writeSource(ofstream& fout, const string& source) const
{
    fout << QTF_RECORD_SOURCE << "\t" << source << endl;
}



void QTFWriter::writeVerticesInQuasiTriangulation(ofstream& fout) const
{
    rg_dList<QTVertex>* q_vertices = m_QT->getVertices();

    QTVertex* currVertex = rg_NULL;
    q_vertices->reset4Loop();
    while ( q_vertices->setNext4Loop() )  {
        currVertex = q_vertices->getpEntity();
        
        if ( currVertex->isVirtual() ) {
            continue;
        }

        rg_Point3D center = currVertex->getBall().getCenter();
        rg_REAL    radius = currVertex->getBall().getRadius();


        fout << QTF_RECORD_VERTEX            << "\t";
        fout << getVertexID( currVertex )  << "\t"; 
        fout << getCellID( currVertex->getFirstTetrahedron() ) << "\t";
        fout.precision(15);
        fout << center.getX()              << "\t";
        fout << center.getY()              << "\t";
        fout << center.getZ()              << "\t";         
        fout << radius                     << "\t";     
        fout << currVertex->getIDFromInput() << endl;
    }
}



void QTFWriter::writeCellsInQuasiTriangulation(ofstream& fout) const
{
    rg_dList<QTTetrahedron>* q_cells = m_QT->getTetrahedra();
    
    q_cells->reset4Loop();
    while ( q_cells->setNext4Loop() )  {
        QTTetrahedron* currCell = q_cells->getpEntity();

        if ( currCell->isVirtual() ) {
            continue;
        }


        QTVertex**      vertices  = currCell->getVertices();
        QTTetrahedron** neighbors = currCell->getNeighbors();

        fout << QTF_RECORD_CELL         << "\t";
        fout << getCellID( currCell ) << "\t";

        rg_INT i;
        for (i=0; i<4; i++) {
            fout << getVertexID( vertices[i] ) << "\t";
        }


        for (i=0; i<4; i++) {
            if ( neighbors[i] != rg_NULL ) {
                fout << getCellID( neighbors[i] ) << "\t";
            }
            else {
                // for isolated triangle.
                fout << "-1" << "\t";
            }
        }

        fout << endl;
    }




    q_cells->reset4Loop();
    while ( q_cells->setNext4Loop() )  {
        QTTetrahedron* currCell = q_cells->getpEntity();

        if ( !currCell->isVirtual() ) {
            continue;
        }


        QTVertex**      vertices  = currCell->getVertices();
        QTTetrahedron** neighbors = currCell->getNeighbors();

        fout << QTF_RECORD_CELL         << "\t";
        fout << getCellID( currCell ) << "\t";

        rg_INT i;
        for (i=0; i<4; i++) {
            fout << getVertexID( vertices[i] ) << "\t";
        }


        for (i=0; i<4; i++) {
            if ( neighbors[i] != rg_NULL ) {
                fout << getCellID( neighbors[i] ) << "\t";
            }
            else {
                // for isolated triangle.
                fout << "-1" << "\t";
            }
        }

        fout << endl;
    }
}



void QTFWriter::writeGatesInQuasiTriangulation(ofstream& fout) const
{
    rg_dList<QTGate>* q_gates = m_QT->getGates()->getGates();


    rg_INT gateID = 1;
    QTGate* currGate = rg_NULL;
    q_gates->reset4Loop();
    while ( q_gates->setNext4Loop() )  {
        currGate = q_gates->getpEntity();

        fout << QTF_RECORD_GATE << "\t";
        fout << gateID        << "\t";
        fout << getVertexID( currGate->getVertex(0) ) << "\t"; 
        fout << getVertexID( currGate->getVertex(1) ) << "\t";


        //  incident big world.
        fout << getCellID( currGate->getBigTetrahedron() ) << "\t";


        //  incident small worlds.
        QTTetrahedron* currWorld = rg_NULL;
        rg_dList<QTTetrahedron*>* smallWorlds = currGate->getSmallTetrahedra();
        smallWorlds->reset4Loop();
        while ( smallWorlds->setNext4Loop() )  {
            currWorld = smallWorlds->getEntity();

            fout << getCellID( currWorld ) << "\t";
        }
        fout << endl;

        gateID++;
    }
}
*/

    

void QTFWriter::writeHeaderForMolecule(ofstream& fout, const string& source, const string& radiusDefinition) 
{
    time_t ltime;
    time( &ltime );
    string creatingTime( ctime( &ltime ) );

    fout << QTF_RECORD_SOURCE        << "\t" << "\"" << source << "\"" << endl;

    if ( !radiusDefinition.empty() ) {
        fout << QTF_RECORD_DEFINE_RADIUS << "\t" << radiusDefinition << endl;
    }

    fout << QTF_RECORD_VDGENENERATOR   << "\t" << "VDRC_VDBALL3D_EDGETRACE_1.0" << endl;
    fout << QTF_RECORD_QTGENENERATOR   << "\t" << "VDRC_1.0" << endl;
    fout << QTF_RECORD_QTFORMATVERSION << "\t" << "1.0" << endl;
    fout << QTF_RECORD_AUTHOR          << "\t" << "VDRC" << "\t" << creatingTime;

    if ( !radiusDefinition.empty() ) {
        fout << QTF_RECORD_COMMENT       << "\t" << "HOH molecules are removed." << endl;
    }
}



void QTFWriter::writeHeaderForPDBMolecule(ofstream& fout, const Molecule* molecule)
{
    time_t ltime;
    time( &ltime );
    string creatingTime( ctime( &ltime ) );

    fout << QTF_RECORD_SOURCE          << "\t" << "\"" << molecule->getMoleculeFileName()    << "\"" << endl;
    fout << QTF_RECORD_PDBTIMESTAMP    << "\t" << molecule->getTimeStamp()                   << endl;
    fout << QTF_RECORD_PDBFILESIZE     << "\t" << molecule->getFileSize()                    << endl;
    fout << QTF_RECORD_DEFINE_RADIUS   << "\t" << molecule->getDescriptionOfAtomRadiusType() << endl;

    fout << QTF_RECORD_VDGENENERATOR   << "\t" << "VDRC_VDBALL3D_EDGETRACE_1.0" << endl;
    fout << QTF_RECORD_QTGENENERATOR   << "\t" << "VDRC_1.0" << endl;
    fout << QTF_RECORD_QTFORMATVERSION << "\t" << "1.0" << endl;
    fout << QTF_RECORD_AUTHOR          << "\t" << "VDRC" << "\t" << creatingTime;

    //fout << QTF_RECORD_COMMENT         << "\t" << "HOH molecules are removed." << endl;
}


   
void QTFWriter::writeQuasiTriangulation(ofstream& fout, const rg_INT& modelID, QuasiTriangulation* QT) 
{
    rg_INT numVertices = QT->getNumBalls();
    rg_INT numCells    = QT->getNumFiniteCells();
    rg_INT numGates    = QT->getNumQTGates();

    fout << QTF_RECORD_MODEL << "\t" 
         << modelID          << "\t"
         << numVertices      << "\t"
         << numCells         << "\t"
         << numGates         << endl;

    writeVerticesInQuasiTriangulation(fout, QT);
    writeCellsInQuasiTriangulation(   fout, QT);
    writeGatesInQuasiTriangulation(   fout, QT);
    writeVDVertexCoordinatesInQuasiTriangulation( fout, QT);

    fout << QTF_RECORD_ENDMODEL << endl;
}



void QTFWriter::writeVerticesInQuasiTriangulation(ofstream& fout, QuasiTriangulation* QT) 
{
    rg_dList<QTVertex>* q_vertices = QT->getVertices();

    QTVertex* currVertex = rg_NULL;
    q_vertices->reset4Loop();
    while ( q_vertices->setNext4Loop() )  {
        currVertex = q_vertices->getpEntity();
        
        if ( currVertex->isVirtual() ) {
            continue;
        }

        rg_Point3D center = currVertex->getBall().getCenter();
        rg_REAL    radius = currVertex->getBall().getRadius();


        fout << QTF_RECORD_VERTEX            << "\t";
        fout << getVertexID( QT, currVertex )  << "\t"; 
        fout << getCellID( currVertex->getFirstTetrahedron() ) << "\t";
        //fout.precision( 3 );
        //fout << fixed << center.getX()              << "\t";
        //fout << fixed << center.getY()              << "\t";
        //fout << fixed << center.getZ()              << "\t";         
        //fout << fixed << radius                     << "\t";     
        fout << center.getX()              << "\t";
        fout << center.getY()              << "\t";
        fout << center.getZ()              << "\t";         
        fout << radius                     << "\t";     
        fout << currVertex->getIDFromInput() << endl;
    }
}



void QTFWriter::writeCellsInQuasiTriangulation(ofstream& fout, QuasiTriangulation* QT) 
{
    rg_dList<QTTetrahedron>* q_cells = QT->getTetrahedra();
    
    q_cells->reset4Loop();
    while ( q_cells->setNext4Loop() )  {
        QTTetrahedron* currCell = q_cells->getpEntity();

        if ( currCell->isVirtual() ) {
            continue;
        }


        QTVertex**      vertices  = currCell->getVertices();
        QTTetrahedron** neighbors = currCell->getNeighbors();

        fout << QTF_RECORD_CELL         << "\t";
        fout << getCellID( currCell ) << "\t";

        rg_INT i;
        for (i=0; i<4; i++) {
            fout << getVertexID( QT, vertices[i] ) << "\t";
        }


        for (i=0; i<4; i++) {
            if ( neighbors[i] != rg_NULL ) {
                if ( neighbors[i]->isVirtual() ) {
                    fout << "0" << "\t";
                }
                else {
                    fout << getCellID( neighbors[i] ) << "\t";
                }
            }
            else {
                // for isolated triangle.
                fout << "-1" << "\t";
            }
        }
        fout << endl;
    }



    /*
    //  write q-cells incident to the virtual vertex.
    q_cells->reset4Loop();
    while ( q_cells->setNext4Loop() )  {
        QTTetrahedron* currCell = q_cells->getpEntity();

        if ( !currCell->isVirtual() ) {
            continue;
        }


        QTVertex**      vertices  = currCell->getVertices();
        QTTetrahedron** neighbors = currCell->getNeighbors();

        fout << QTF_RECORD_CELL         << "\t";
        fout << getCellID( currCell ) << "\t";

        rg_INT i;
        for (i=0; i<4; i++) {
            fout << getVertexID( QT, vertices[i] ) << "\t";
        }


        for (i=0; i<4; i++) {
            if ( neighbors[i] != rg_NULL ) {
                fout << getCellID( neighbors[i] ) << "\t";
            }
            else {
                // for isolated triangle.
                fout << "-1" << "\t";
            }
        }

        fout << endl;
    }
    */
}



void QTFWriter::writeGatesInQuasiTriangulation(ofstream& fout, QuasiTriangulation* QT) 
{
    rg_dList<QTGate>* q_gates = QT->getGates()->getGates();


    rg_INT gateID = 1;
    QTGate* currGate = rg_NULL;
    q_gates->reset4Loop();
    while ( q_gates->setNext4Loop() )  {
        currGate = q_gates->getpEntity();

        fout << QTF_RECORD_GATE << "\t";
        fout << gateID        << "\t";
        fout << getVertexID( QT, currGate->getVertex(0) ) << "\t"; 
        fout << getVertexID( QT, currGate->getVertex(1) ) << "\t";


        //  incident big world.
        fout << getCellID( currGate->getBigTetrahedron() ) << "\t";


        //  incident small worlds.
        QTTetrahedron* currWorld = rg_NULL;
        rg_dList<QTTetrahedron*>* smallWorlds = currGate->getSmallTetrahedra();
        smallWorlds->reset4Loop();
        while ( smallWorlds->setNext4Loop() )  {
            currWorld = smallWorlds->getEntity();

            fout << getCellID( currWorld ) << "\t";
        }
        fout << endl;

        gateID++;
    }
}



void QTFWriter::writeVDVertexCoordinatesInQuasiTriangulation(   ofstream& fout, QuasiTriangulation* QT)
{
    rg_dList<QTTetrahedron>* q_cells = QT->getTetrahedra();
    
    q_cells->reset4Loop();
    while ( q_cells->setNext4Loop() )  {
        QTTetrahedron* currCell = q_cells->getpEntity();

        if ( currCell->isVirtual() ) {
            continue;
        }

        if ( currCell->isIsolatedFace() ) {
            continue;
        }


        fout << QTF_RECORD_VDVTX      << "\t";
        fout << getCellID( currCell ) << "\t";

        rg_Point3D vdVtxCoord = currCell->getEmptyTangentSphere().getCenter();

        fout.precision( 4 );
        //fout << vdVtxCoord.getX() << "\t";
        //fout << vdVtxCoord.getY() << "\t";
        //fout << vdVtxCoord.getZ() << endl;
        fout << fixed << vdVtxCoord.getX() << "\t";
        fout << fixed << vdVtxCoord.getY() << "\t";
        fout << fixed << vdVtxCoord.getZ() << endl;
    }
}



