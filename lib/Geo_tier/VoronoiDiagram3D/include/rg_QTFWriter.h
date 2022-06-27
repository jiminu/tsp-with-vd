#ifndef _QTFWRITER_H
#define _QTFWRITER_H

#include "rg_Const.h"
#include "rg_QuasiTriangulation.h"
#include "rg_Molecule.h"
#include <fstream>
#include <string>
using namespace std;

namespace V {

namespace GeometryTier {


class QTFWriter
{
private:
    string  m_QTFileName;
    QuasiTriangulation* m_QT;

public:
    QTFWriter();
    QTFWriter(QuasiTriangulation* QT);
    QTFWriter(const QTFWriter& writer);
    ~QTFWriter();

    QuasiTriangulation* getQuasiTriangulation() const;
    void                setQuasiTriangulation(QuasiTriangulation* QT);

    QTFWriter&      operator=(const QTFWriter& writer);


    void write( const string& qtfFilename, const string& source, QuasiTriangulation* QT, Molecule* molecule );
    
    void write( const string& qtfFilename, 
                const string& source,
                const rg_dList< pair<QuasiTriangulation*, Molecule*> >& models );
    void write( const string& qtfFilename, const string& source = "");

private:
    rg_INT getVertexID(QuasiTriangulation* QT, QTVertex* vertex);
    rg_INT getCellID(QTTetrahedron* cell);

    void writeHeaderForMolecule(ofstream& fout, const string& source, const string& radiusDefinition);
    void writeHeaderForPDBMolecule(ofstream& fout, const Molecule* molecule);

    void writeQuasiTriangulation(ofstream& fout, const rg_INT& modelID, QuasiTriangulation* QT);
        void writeVerticesInQuasiTriangulation(ofstream& fout, QuasiTriangulation*  QT);
        void writeCellsInQuasiTriangulation(   ofstream& fout, QuasiTriangulation* QT);
        void writeGatesInQuasiTriangulation(   ofstream& fout, QuasiTriangulation* QT);
        void writeVDVertexCoordinatesInQuasiTriangulation(   ofstream& fout, QuasiTriangulation* QT);

//     void writeSource(ofstream& fout, const string& source) const;
//     void writeVerticesInQuasiTriangulation(ofstream& fout) const;
//     void writeCellsInQuasiTriangulation(ofstream& fout) const;
//     void writeGatesInQuasiTriangulation(ofstream& fout) const;


};

} // namespace GeometryTier

} // namespace V


#endif
