#ifndef _VDCELL_H
#define _VDCELL_H

#include "ConstForVoronoiDiagram3D.h"

#include "rg_dList.h"

#include "TopologicalEntity.h"

//#include <fstream>
//using namespace std;

// aaa

namespace V {

namespace GeometryTier {


class BallGenerator;
class VDEdge;
class VDFace;
class VDVertex;

class QTVertex;
class BetaVertex;


class VDCell : public TopologicalEntity
{
private:
    rg_dList<VDFace*> m_boundingFaces;

    BallGenerator*    m_generator;

    rg_FLAG m_boundedness;

    QTVertex*   m_qVertex;
    BetaVertex* m_betaVertex;

public:
    //  constructor & deconstructor..
    VDCell();
    VDCell(const rg_INT& ID);
    VDCell(BallGenerator* generator);
    VDCell(const rg_INT& ID, BallGenerator* generator);
    VDCell(const TopologicalEntity& aTopoEntity, BallGenerator* generator);
    VDCell(const VDCell& aCell);
    ~VDCell();

    //  get functions.. 
    inline QTVertex* getQVertex() const { return m_qVertex; }
    inline BetaVertex* getBetaVertex() const { return m_betaVertex; }


    BallGenerator*     getGenerator() const;
    rg_dList<VDFace*>* getBoungindFaces();
    rg_INT             getNumOfBoundingFaces() const;

    rg_FLAG isBounded() const;
    void    isBounded( const rg_FLAG& boundedness);

    //  set functions..
    void setGenerator(BallGenerator* generator);
    void setBoundingFaces(const rg_dList<VDFace*>& boundingFaces);
    void addBoundingFace(VDFace* aBoundingFace);

    friend void connectCellAndGenerator(VDCell* cell, BallGenerator* generator);

	//computation
	rg_REAL computeVolume();
	rg_REAL computeArea();

    //  operator overloading..
    VDCell& operator =(const VDCell& aCell);

    //  topological operators..
    VDFace* findFaceToShareWith(VDCell* anotherCell);
    rg_INT  findFacesToShareWith(VDCell* anotherCell, rg_dList<VDFace*>& faceList) const;

    //  topological quires..
    void  inquireBoundingVertices(rg_dList<VDVertex*>& vertexList);
    void  inquireBoundingEdges(rg_dList<VDEdge*>& edgeList);
    void  inquireBoundingFaces(rg_dList<VDFace*>& faceList);
    void  inquireIncidentCells(rg_dList<VDCell*>& cellList);


    //  connect vertices in IWDS and eIWDS
    void connectQVertex(   QTVertex*   q_vtx);
    void disconnectQVertex(QTVertex*   q_vtx);

    void connectBetaVertex(   BetaVertex* b_vtx);
    void disconnectBetaVertex(BetaVertex* b_vtx);
};


void connectCellAndGenerator(VDCell* cell, BallGenerator* generator);


} // namespace GeometryTier

} // namespace V


#endif

