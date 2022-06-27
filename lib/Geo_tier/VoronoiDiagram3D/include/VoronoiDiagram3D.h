#ifndef _VORONOIDIAGRAM3D_H
#define _VORONOIDIAGRAM3D_H

#include <fstream>
using namespace std;

#include "ConstForVoronoiDiagram3D.h"

#include "rg_dList.h"

#include "VDCell.h"
#include "VDFace.h"
#include "VDLoop.h"
#include "VDPartialEdge.h"
#include "VDEdge.h"
#include "VDVertex.h"

#include "Bucket.h"
#include <map>
using namespace std;

namespace V {

namespace GeometryTier {


class VoronoiDiagram3D
{
protected:

    rg_dList< VDCell >        m_GlobalCellList;
    rg_dList< VDFace >        m_GlobalFaceList;
    rg_dList< VDLoop >        m_GlobalLoopList;
    rg_dList< VDPartialEdge > m_GlobalPartialedgeList;
    rg_dList< VDEdge >        m_GlobalEdgeList;
    rg_dList< VDVertex >      m_GlobalVertexList;


    Bucket m_bucket;
//    Bucket m_bucket2;

    rg_INT     m_status;

public:
    rg_FLAG    m_makeVoronoiFace;

    //  constructor & deconstructor..
    VoronoiDiagram3D();
    ~VoronoiDiagram3D();

    void duplicateTopology(const VoronoiDiagram3D& origin);
    void clean();
    void cleanForRestart();

    //  get functions.. 
    rg_dList< VDCell >*        getGlobalCellList();
    rg_dList< VDFace >*        getGlobalFaceList();
    rg_dList< VDLoop >*        getGlobalLoopList();
    rg_dList< VDPartialEdge >* getGlobalPartialEdgeList();
    rg_dList< VDEdge >*        getGlobalEdgeList();
    rg_dList< VDVertex >*      getGlobalVerticeList();

    inline const rg_dList< VDCell >&        getGlobalCellList() const        { return m_GlobalCellList; }
    inline const rg_dList< VDFace >&        getGlobalFaceList() const        { return m_GlobalFaceList; }
    inline const rg_dList< VDLoop >&        getGlobalLoopList() const        { return m_GlobalLoopList; }
    inline const rg_dList< VDPartialEdge >& getGlobalPartialEdgeList() const { return m_GlobalPartialedgeList; }
    inline const rg_dList< VDEdge >&        getGlobalEdgeList() const        { return m_GlobalEdgeList; }
    inline const rg_dList< VDVertex >&      getGlobalVerticeList() const     { return m_GlobalVertexList; }

    rg_INT getNumOfCells() const;
    rg_INT getNumOfFaces() const;
    rg_INT getNumOfLoops() const;
    rg_INT getNumOfPartialEdges() const;
    rg_INT getNumOfEdges() const;
    rg_INT getNumOfVertices() const;

    rg_INT getNumOfDisconnectedGenerator();

    rg_INT getStatus() const;
    
    //  set functions..



    //  operator overloading..
    void fileOutGlobalEdgeList( ofstream& fout );
    void fileOutEdgeLengthsOfDelaunayTetrahedron( ofstream& fout );

    void writeCellList( ofstream& fout );
    void writeFaceList( ofstream& fout );
    void writeLoopList( ofstream& fout );
    void writePrEdgeList( ofstream& fout );
    void writeEdgeList( ofstream& fout );
    void writeVertexList( ofstream& fout );
    void writeGeneratorList( ofstream& fout );

    //  topological operators..
    VDVertex*       createVertex(const rg_FLAG& onInfinity, const rg_Point3D& point, const rg_REAL& radiusOfTS);
    VDEdge*         createEdge(VDVertex* startVtx, VDVertex* endVtx);
    VDPartialEdge*  createPrEdge(VDEdge* oriEdge);


    VDCell*        createCell(  const VDCell& cell);
    VDFace*        createFace(  const VDFace& face);
    VDLoop*        createLoop(  const VDLoop& loop);
    VDPartialEdge* createPrEdge(const VDPartialEdge& prEdge);
    VDEdge*        createEdge(  const VDEdge& edge);
    VDVertex*      createVertex(const VDVertex& vertex);


    void splitEdge(VDEdge* edge, VDVertex* splitVtx, VDEdge*& newEdge1, VDEdge*& newEdge2);


    VDVertex* findClosestVertexToGivenPoint( const rg_Point3D& givenPoint );


    //  geometric operators..
    rg_dList<BallGenerator*>* computeConvexHullOfGeneter();
    rg_dList<VDFace*>*    obtainFacesWithIntersectingGenerators();
    rg_dList<VDEdge*>*    obtainEdgesWithoutIntersectingGenerators();
    rg_dList<VDCell*>*    obtainCellsToIntersectWithSectionalYZPlane(const rg_REAL& xOfPlane);
    
    void computeVoronoiFaces();
    void computeVoronoiFaces(rg_INT& numTriangles);

    rg_INT approximateStorageCost();


private:
    void makeMapForDuplication( const VoronoiDiagram3D& origin, 
                                map<VDCell*, VDCell*>&               cellMap,
                                map<VDFace*, VDFace*>&               faceMap,
                                map<VDLoop*, VDLoop*>&               loopMap,
                                map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                map<VDEdge*, VDEdge*>&               edgeMap,
                                map<VDVertex*, VDVertex*>&           vertexMap);
    void duplicateTopologyOfCells( const VoronoiDiagram3D& origin, 
                                   map<VDCell*, VDCell*>&               cellMap,
                                   map<VDFace*, VDFace*>&               faceMap);
    void duplicateTopologyOfFaces( const VoronoiDiagram3D& origin, 
                                   map<VDCell*, VDCell*>&               cellMap,
                                   map<VDFace*, VDFace*>&               faceMap,
                                   map<VDLoop*, VDLoop*>&               loopMap);
    void duplicateTopologyOfLoops( const VoronoiDiagram3D& origin, 
                                   map<VDFace*, VDFace*>&               faceMap,
                                   map<VDLoop*, VDLoop*>&               loopMap,
                                   map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap);
    void duplicateTopologyOfPartialEdges( const VoronoiDiagram3D& origin, 
                                   map<VDLoop*, VDLoop*>&               loopMap,
                                   map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                   map<VDEdge*, VDEdge*>&               edgeMap);
    void duplicateTopologyOfEdges( const VoronoiDiagram3D& origin, 
                                   map<VDPartialEdge*, VDPartialEdge*>& prEdgeMap,
                                   map<VDEdge*, VDEdge*>&               edgeMap,
                                   map<VDVertex*, VDVertex*>&           vertexMap);
    void duplicateTopologyOfVertices( const VoronoiDiagram3D& origin, 
                                   map<VDEdge*, VDEdge*>&               edgeMap,
                                   map<VDVertex*, VDVertex*>&           vertexMap);

};

} // namespace GeometryTier

} // namespace V


#endif
