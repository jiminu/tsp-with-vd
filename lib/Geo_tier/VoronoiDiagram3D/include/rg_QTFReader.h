#ifndef _QTFREADER_H
#define _QTFREADER_H

#include "rg_Const.h"
#include "rg_dList.h"
#include "rg_QuasiTriangulation.h"
#include "rg_Molecule.h"

#include <string>
#include <list>
#include <fstream>
using namespace std;



namespace V {

namespace GeometryTier {


class QTFReader
{
private:
    QuasiTriangulation* m_QT;


    struct QTIncexVertex {
        QTVertex*      qtVertex;
        rg_INT         firstCell;
    };

    struct QTIndexCell {
        QTTetrahedron* qtCell;
        rg_INT         vertexIndex[4];
        rg_INT         neighborIndex[4];
    };


    rg_BOOL                   m_useVDVerticesInQTFile;
    rg_dList< QTIncexVertex > m_indexVertices;
    rg_dList< QTIndexCell >   m_indexCells;

    QTVertex**      m_vertices;
    QTTetrahedron** m_cells;


public:
    //// for analysis of QTFReader
    //double timeLoading;
    //double timeParsing;
    //double timeSettingQT;
    //double timeParsingVVtx;
    //double timeComputingVVtx;

    QTFReader();
    QTFReader(QuasiTriangulation* QT);
    QTFReader(const QTFReader& reader);
    ~QTFReader();

    QuasiTriangulation* getQuasiTriangulation() const;
    void                setQuasiTriangulation(QuasiTriangulation* QT);

    QTFReader&      operator=(const QTFReader& reader);

    rg_BOOL read(const Molecule* const& molecule, const string& qtfFilename, QuasiTriangulation* m_QT, const rg_BOOL& useVDVertices = rg_TRUE);
    rg_BOOL read(const string& qtfFilename, const rg_INT& modelID = 1, const rg_BOOL& useVDVertices = rg_TRUE);
    rg_BOOL read2(const string& qtfFilename, const rg_INT& modelID = 1, const rg_BOOL& useVDVertices = rg_TRUE);

    void computeEmptyTangentSpheresOfCells();

private:
    string parsePDBTimeStamp(const string& qtfFilename);
    rg_INT parsePDBFileSize(const string& qtfFilename);


    void extractSingleModelFromQTFile(const string& qtfFilename, const rg_INT& modelID, list<string>& recordsOfSelectedModel );

    void parseVerticesOfQTFile(list<string>&  recordsInModel,  list<string*>& recordsForCells );
    void parseCellsOfQTFile(   list<string*>& recordsForCells, list<string*>& recordsForGates );

    void setTopologyBetweenVerticesAndCells();

    void parseGatesOfQTFile(   list<string*>& recordsForGates, list<string*>& recordsForVDVertices);
    void parseVDVerticesOfQTFile(list<string*>& recordsForVDVertices);


    void extractSingleModelFromQTFile( const string& qtfFilename, const rg_INT& modelID, 
                                              list<string>& recordsForQTVertices,
                                              list<string>& recordsForQTCells,
                                              list<string>& recordsForQTGates,
                                              list<string>& recordsForVDVertices );
    void parseVerticesOfQTFile(const list<string>& recordsForQTVertices );
    void parseCellsOfQTFile(   const list<string>& recordsForQTCells );
    void parseGatesOfQTFile(   const list<string>& recordsForQTGates );
    void parseVDVerticesOfQTFile(const list<string>& recordsForVDVertices);

    rg_BOOL synchronizeAtomsBetweenQTAndMolecule( const Molecule* const& molecule );


    //void extractQuasiTriangulation(list<string>& recordsOfQTFFile, const rg_INT& modelID, list<string*>& recordsOfQTInModel);
    //void readVerticesInQuasiTriangulation(ifstream& fin, list<string>& recordsOfQTFFile);
    //void readCellsInQuasiTriangulation(ifstream& fin, list<string>& recordsOfQTFFile);
    //void readGatesInQuasiTriangulation(ifstream& fin, list<string>& recordsOfQTFFile);
};

} // namespace GeometryTier

} // namespace V


#endif
