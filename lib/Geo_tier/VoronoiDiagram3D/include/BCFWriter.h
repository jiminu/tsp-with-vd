#ifndef BCFWRITER_H
#define BCFWRITER_H

#include "rg_Const.h"
#include "BetaUniverse.h"
#include "rg_Molecule.h"

#include <map>
using namespace std;



namespace V {

namespace GeometryTier {


class BCFWriter
{
public:
    BCFWriter();
    ~BCFWriter();

    static void write(const string& bcfFilename, const string& source, Molecule* molecule, BetaUniverse* QT, const rg_REAL& betaValue, const rg_BOOL& onlyBetaShape=rg_FALSE);

private:
    static void writeHeaderForPDBMolecule(ofstream& fout, const Molecule* molecule, const rg_REAL& betaValue, const rg_BOOL& onlyBetaShape);
    static void writeHeaderForMolecule(ofstream& fout, const string& source, const string& radiusDefinition);
    static void writeBetaComplex(ofstream& fout, const rg_INT& modelID, BetaUniverse* QT, const rg_REAL& betaValue);
    static void writeBetaShape(ofstream& fout, const rg_INT& modelID, BetaUniverse* QT, const rg_REAL& betaValue);
        static void writeVertices(ofstream& fout, const rg_dList<BetaVertex*>& vertices, const rg_REAL& betaValue, map<BetaVertex*, rg_INT>& vertexIDs);
        static void writeEdges(ofstream& fout, const rg_dList<BetaEdge*>& edges, const rg_REAL& betaValue, const map<BetaVertex*, rg_INT>& vertexIDs);
        static void writeFaces(ofstream& fout, const rg_dList<BetaFace*>& faces, const rg_REAL& betaValue, const map<BetaVertex*, rg_INT>& vertexIDs);
        static void writeCells(ofstream& fout, const rg_dList<BetaCell*>& cells, const rg_REAL& betaValue, const map<BetaVertex*, rg_INT>& vertexIDs);

    static void writeStatisticsOfBetaComplexOrShape(ofstream& fout, BetaUniverse* QT, const rg_REAL& betaValue, const rg_BOOL& onlyBetaShape);

    inline static bool BCVertexInputIDLess(BetaVertex* vertex1, BetaVertex* vertex2) {
        return vertex1->getBallProperty()->getIDFromInput() < vertex2->getBallProperty()->getIDFromInput();
    }

public:
    class BetaEdgeByVertexID
    {
    private:
        BetaEdge* m_edge;
        rg_INT    m_vtxIDs[2];
    public:
        BetaEdgeByVertexID();
        BetaEdgeByVertexID(BetaEdge* edge, const rg_INT& vtxID1, const rg_INT& vtxID2);
        BetaEdgeByVertexID(const BetaEdgeByVertexID& betaEdge);
        ~BetaEdgeByVertexID();

        inline BetaEdge* getEdge() const { return m_edge; }
        inline rg_INT    getVertexID(const rg_INT& i) const { return m_vtxIDs[i]; }
        void      setBetaEdgeByVertexID(BetaEdge* edge, const rg_INT& vtxID1, const rg_INT& vtxID2);
        BetaEdgeByVertexID& operator =(const BetaEdgeByVertexID& betaEdge);
        bool operator==(const BetaEdgeByVertexID& betaEdge) const;
        bool operator <(const BetaEdgeByVertexID& betaEdge) const;
    };



    class BetaFaceByVertexID
    {
    private:
        BetaFace* m_face;
        rg_INT    m_vtxIDs[3];
    public:
        BetaFaceByVertexID();
        BetaFaceByVertexID(BetaFace* face, const rg_INT& vtxID1, const rg_INT& vtxID2, const rg_INT& vtxID3);
        BetaFaceByVertexID(const BetaFaceByVertexID& betaFace);
        ~BetaFaceByVertexID();

        inline BetaFace* getFace() const { return m_face; }
        inline rg_INT    getVertexID(const rg_INT& i) const { return m_vtxIDs[i]; }
        void      setBetaFaceByVertexID(BetaFace* face, const rg_INT& vtxID1, const rg_INT& vtxID2, const rg_INT& vtxID3);
        BetaFaceByVertexID& operator =(const BetaFaceByVertexID& betaFace);
        bool operator==(const BetaFaceByVertexID& betaFace) const;
        bool operator <(const BetaFaceByVertexID& betaFace) const;
    };


    class BetaCellByVertexID
    {
    private:
        BetaCell* m_cell;
        rg_INT    m_vtxIDs[4];
        bool      m_positiveVolume;
    public:
        BetaCellByVertexID();
        BetaCellByVertexID(BetaCell* cell, const rg_INT& vtxID1, const rg_INT& vtxID2, const rg_INT& vtxID3, const rg_INT& vtxID4, const bool& positiveVolume);
        BetaCellByVertexID(const BetaCellByVertexID& betaCell);
        ~BetaCellByVertexID();

        inline BetaCell* getFace() const { return m_cell; }
        inline rg_INT    getVertexID(const rg_INT& i) const { return m_vtxIDs[i]; }
        inline bool      hasPositiveVolume() const { return m_positiveVolume; }
        void      setBetaCellByVertexID(BetaCell* cell, const rg_INT& vtxID1, const rg_INT& vtxID2, const rg_INT& vtxID3, const rg_INT& vtxID4, const bool& positiveVolume);
        BetaCellByVertexID& operator =(const BetaCellByVertexID& betaCell);
        bool operator==(const BetaCellByVertexID& betaCell) const;
        bool operator <(const BetaCellByVertexID& betaCell) const;
    };

};


} // namespace GeometryTier

} // namespace V




#endif 


