#ifndef _VDFACE_H
#define _VDFACE_H

#include "ConstForVoronoiDiagram3D.h"

#include "rg_dList.h"
#include "rg_Surface.h"

#include "TopologicalEntity.h"
#include "Sphere.h"

//#include "tri_header.h"

#include "TriangleOnVDFace.h"

#include <fstream>
using namespace std;

namespace V {

namespace GeometryTier {


class VDCell;
class VDLoop;
class VDEdge;
class VDVertex;

class BetaEdge;

const int numOfEvalPtsOnCrv = 5;
const int numOfGridPts = 5;

const int NUM_POINTS_ON_CURVE = 20;
const int NUM_GRID_POINTS     = 50;

class VDFace : public TopologicalEntity
{
private:
    VDCell* m_rightCell;
    VDCell* m_leftCell;

    rg_dList<VDLoop*> m_loops;


    rg_Surface*           m_faceEquation;
    rg_dList<rg_Point3D>* m_faceMesh;

    rg_dList<TriangleOnVDFace>* m_triangularMesh;


    rg_FLAG    m_isOnInfinity;
    rg_FLAG    m_ProbeTangibility;
    rg_FLAG    m_boundedness;


    rg_BOOL    m_visited;
    BetaEdge*  m_betaEdge;


public:
    //  constructor & deconstructor..
    VDFace();
    VDFace(const rg_INT& ID); 
    VDFace(const rg_INT& ID, VDCell* rightCell, VDCell* leftCell);
    VDFace(const TopologicalEntity& aTopoEntity, VDCell* rightCell, VDCell* leftCell);
    VDFace(const VDFace& aFace);
    ~VDFace();

    //  get functions.. 
    inline BetaEdge* getBetaEdge() const { return m_betaEdge; }
    inline rg_BOOL   isVisited()   const { return m_visited; }
    inline void      isVisited(const rg_BOOL& visited) { m_visited = visited; }

    VDCell* getRightCell() const;
    VDCell* getLeftCell() const;
    
    VDLoop* getOuterLoop() const;
    VDLoop* getLoop(const rg_INT& i) const;
    rg_INT  getNumOfLoops() const;

    rg_dList<VDLoop*>* getLoops();

    rg_Surface*           getFaceEquation() const;
    rg_dList<rg_Point3D>* getFaceMesh();

    rg_dList<TriangleOnVDFace>* getTriangularMesh();

    rg_BOOL isVirtual() const;
    rg_BOOL isUnbounded() const;


    rg_FLAG isOnInfinity() const;
	void setInfinity();

    rg_FLAG isTangible() const;
    void    isTangible(const rg_FLAG& tangibility);

    rg_FLAG isBounded() const;
    void    isBounded( const rg_FLAG& boundedness);

    rg_FLAG isThereInnerLoop();

    //  set functions..
    void setRightCell(VDCell* rightCell);
    void setLeftCell(VDCell* leftCell);

    void connectFaceWithRightCell(VDCell* rightCell);
    void connectFaceWithLeftCell(VDCell* leftCell);
    void connectFaceWithRightAndLeftCell(VDCell* rightCell, VDCell* leftCell);


    void setLoops(const rg_dList<VDLoop*>& loops);
    void addLoop(VDLoop* loop);


    void setFaceEquation(rg_Surface* faceEquation);

    //  operator overloading..
	VDFace& operator=(const VDFace& aFace);

    /*

	//computation
	rg_REAL computeArea();

    //  topological operators..

    //  geometric operators..
    rg_FLAG isThereIntersectionWithSectionalYZPlane(const rg_REAL& xOfPlane);

    void makeMeshForVoronoiFace();
    void makeTriangularMeshForVoronoiFace();
        rg_REAL getMaxDistanceBetweenTrianglesAndGenerator();

//    rg_dList<rg_Point3D>* triangulateVoronoiFace(const Sphere& sphere1, const Sphere& sphere2, CrvList& bzCrvList);
    rg_dList<rg_Point3D>* triangulateVoronoiFace(const Sphere& sphere1, const Sphere& sphere2, CrvList& bzCrvList, 
                                                rg_INT numPointsOnCurve = NUM_POINTS_ON_CURVE, 
                                                rg_INT numGridPoints    = NUM_GRID_POINTS);
//    void makeMeshForVoronoiFace(ofstream& fout);
//    rg_dList<rg_Point3D>* triangulateVoronoiFace(const Sphere& sphere1, const Sphere& sphere2, CrvList& bzCrvList, ofstream& fout);
        Point2d    centerOfMass(Face_iterator& it);
        Point3d    getPoint3d(const rg_Point3D& pt);
        rg_REAL    zValue(const leda_real& leda_x, const leda_real& leda_y, const rg_REAL& f, const rg_REAL& r);
		rg_REAL    zValue(const rg_REAL& x, const rg_REAL& y, const rg_REAL& f, const rg_REAL& r);
        rg_Point3D getRGPoint3D(const Point3d& pt);
        Point3d    getNormal(const Point3d& pt, const rg_REAL& f,const rg_REAL& r);
    */

    //  topological quires..
    void  inquireBoundingVertices(rg_dList<VDVertex*>& vertexList) const;
    void  inquireBoundingEdges(rg_dList<VDEdge*>& edgeList) const;
    void  inquireIncidentFaces(rg_dList<VDFace*>& faceList);
    void  inquireIncidentCells(rg_dList<VDCell*>& cellList);

    rg_BOOL isIncidentTo(   VDCell* currCell ) const;



    void connectBetaEdge(   BetaEdge* b_edge);
    void disconnectBetaEdge(BetaEdge* b_edge);
};


} // namespace GeometryTier

} // namespace V


#endif

