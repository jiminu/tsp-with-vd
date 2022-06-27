#ifndef _QUASITRIANGULATION_H
#define _QUASITRIANGULATION_H

#include "rg_dList.h"
#include "rg_QTVertex.h"
#include "rg_QTTetrahedron.h"
#include "rg_QTGate.h"
#include "rg_QTGateList.h"

#include "Ball.h"

#include "rg_Circle2D.h"
#include "rg_ImplicitEquation.h"

#include <fstream>
using namespace std;

namespace V {

namespace GeometryTier {



const int DEALWITH_NO_QTENTITY_IN_QTF   = 0;
const int DEALWITH_QTVERTEX_IN_QTF      = 1;
const int DEALWITH_QTTETRAHEDRON_IN_QTF = 2;
const int DEALWITH_QTGATE_IN_QTF        = 3;



class QTEdge;
class QTFace;

class QuasiTriangulation
{
private:
    rg_dList< Ball >          m_balls;

    rg_dList< QTVertex >      m_vertices;
    rg_dList< QTTetrahedron > m_tetrahedra;
    QTGateList                m_gates;

    rg_BOOL                   m_isLinkedVD;

public:
    QuasiTriangulation();
    ~QuasiTriangulation();

    rg_INT getNumBalls() const;
    rg_INT getNumQTVertices() const;
    rg_INT getNumQTTetrahedron() const;
    rg_INT getNumQTGates() const;

    rg_INT getNumFiniteCells();

    rg_dList<QTVertex>*      getVertices();
    rg_dList<QTTetrahedron>* getTetrahedra();
    QTGateList*              getGates();

    rg_INT         getAbsoluteIDOfVertex(QTVertex* vertex);

    Ball*          addBall(const Ball& ball);
    QTVertex*      addVertex(const QTVertex& vertex);
    QTTetrahedron* addTetrahedron(const QTTetrahedron& tetrahedron);
    QTGate*        addGate(const QTGate& gate);


    rg_FLAG isOnConvexHull(QTVertex* vertex);
    rg_FLAG isOnConvexHull(const QTEdge& edge);
    rg_FLAG isOnConvexHull(const QTFace& face);
    rg_FLAG isOnConvexHull(QTTetrahedron* tetrahedron);

    void    resetIDs();

    inline rg_BOOL isLinkedVD() const { return m_isLinkedVD; }
    inline void    isLinkedVD(const rg_BOOL& bLinkedVD) { m_isLinkedVD = bLinkedVD; }

    ///////////////////////////////////////////////////////////////////////////////
    //    For Queres  
    //   T -> V, E, F, T
    void inquireBoundingVerticesOfTetrahedronInThisWorld(QTTetrahedron* givenTetrahedron, rg_dList<QTVertex*>& vertexList);
    void inquireBoundingEdgesOfTetrahedronInThisWorld(QTTetrahedron* givenTetrahedron, rg_dList<QTEdge>& edgeList);
    void inquireBoundingFacesOfTetrahedronInThisWorld(QTTetrahedron* givenTetrahedron, rg_dList<QTFace>& faceList);
    void inquireIncidentTetrahedraOfTetrahedronInThisWorld(QTTetrahedron* givenTetrahedron, rg_dList<QTTetrahedron*>& tetrahedronList);

    //   F -> V, E, F, T
    void inquireBoundingVerticesOfFaceInThisWorld(const QTFace& givenFace, rg_dList<QTVertex*>& vertexList);
    void inquireBoundingEdgesOfFaceInThisWorld(const QTFace& givenFace, rg_dList<QTEdge>& edgeList);
    void inquireIncidentFacesOfFaceInThisWorld(const QTFace& givenFace, rg_dList<QTFace>& faceList);
    void inquireIncidentTetrahedraOfFaceInThisWorld(const QTFace& givenFace, rg_dList<QTTetrahedron*>& tetrahedronList);

    //   E -> V, E, F, T
    void inquireBoundingVerticesOfEdgeInThisWorld(const QTEdge& givenEdge, rg_dList<QTVertex*>& vertexList);
    void inquireIncidentEdgesOfEdgeInThisWorld(const QTEdge& givenEdge, rg_dList<QTEdge>& edgeList);
    void inquireIncidentFacesOfEdgeInThisWorld(const QTEdge& givenEdge, rg_dList<QTFace>& faceList);
    void inquireIncidentTetrahedraOfEdgeInThisWorld(const QTEdge& givenEdge, rg_dList<QTTetrahedron*>& tetrahedronList);

    //   V -> V, E, F, T
    void inquireNeighborVerticesOfVertexInThisWorld(QTVertex* givenVertex, rg_dList<QTVertex*>& vertexList);
    void inquireIncidentEdgesOfVertexInThisWorld(QTVertex* givenVertex, rg_dList<QTEdge>& edgeList);
    void inquireIncidentFacesOfVertexInThisWorld(QTVertex* givenVertex, rg_dList<QTFace>& faceList);
    void inquireIncidentTetrahedraOfVertexInThisWorld(QTVertex* givenVertex, rg_dList<QTTetrahedron*>& tetrahedronList);


    ///////////////////////////////////////////////////////////////////////////
    //
    //  Quasi-operator
    //
    //  1. Query primitives for intra-world.
    //
    //  1.1. Q(v, Y | W): Query primitives for intra-world.
    void inquireVerticesOfVertexInIntraWorld(QTVertex* givenVertex, rg_dList<QTVertex*>& vertexList,           QTTetrahedron* world = rg_NULL);
    void inquireEdgesOfVertexInIntraWorld(   QTVertex* givenVertex, rg_dList<QTEdge>& edgeList,               QTTetrahedron* world = rg_NULL);
    void inquireFacesOfVertexInIntraWorld(   QTVertex* givenVertex, rg_dList<QTFace>& faceList,               QTTetrahedron* world = rg_NULL);
    void inquireCellsOfVertexInIntraWorld(   QTVertex* givenVertex, rg_dList<QTTetrahedron*>& tetrahedronList, QTTetrahedron* world = rg_NULL);
    
    //  1.2. Q(e, Y | W): Query primitives for intra-world.
    void inquireVerticesOfEdgeInIntraWorld(const QTEdge& givenEdge, rg_dList<QTVertex*>& vertexList,           QTTetrahedron* world = rg_NULL);
    void inquireEdgesOfEdgeInIntraWorld(   const QTEdge& givenEdge, rg_dList<QTEdge>& edgeList,               QTTetrahedron* world = rg_NULL);
    void inquireFacesOfEdgeInIntraWorld(   const QTEdge& givenEdge, rg_dList<QTFace>& faceList,               QTTetrahedron* world = rg_NULL);
    void inquireCellsOfEdgeInIntraWorld(   const QTEdge& givenEdge, rg_dList<QTTetrahedron*>& tetrahedronList, QTTetrahedron* world = rg_NULL);
    
    //  1.3. Q(f, Y | W): Query primitives for intra-world.
    void inquireVerticesOfFaceInIntraWorld(const QTFace& givenFace, rg_dList<QTVertex*>& vertexList,           QTTetrahedron* world = rg_NULL);
    void inquireEdgesOfFaceInIntraWorld(   const QTFace& givenFace, rg_dList<QTEdge>& edgeList,               QTTetrahedron* world = rg_NULL);
    void inquireFacesOfFaceInIntraWorld(   const QTFace& givenFace, rg_dList<QTFace>& faceList,               QTTetrahedron* world = rg_NULL);
    void inquireCellsOfFaceInIntraWorld(   const QTFace& givenFace, rg_dList<QTTetrahedron*>& tetrahedronList, QTTetrahedron* world = rg_NULL);
    
    //  1.4. Q(c, Y | W): Query primitives for intra-world.
    void inquireVerticesOfCellInIntraWorld(QTTetrahedron* givenCell, rg_dList<QTVertex*>& vertexList,           QTTetrahedron* world = rg_NULL);
    void inquireEdgesOfCellInIntraWorld(   QTTetrahedron* givenCell, rg_dList<QTEdge>& edgeList,               QTTetrahedron* world = rg_NULL);
    void inquireFacesOfCellInIntraWorld(   QTTetrahedron* givenCell, rg_dList<QTFace>& faceList,               QTTetrahedron* world = rg_NULL);
    void inquireCellsOfCellInIntraWorld(   QTTetrahedron* givenCell, rg_dList<QTTetrahedron*>& tetrahedronList, QTTetrahedron* world = rg_NULL);
    
    //  2. Query primitives for inter-world.
    //
    void inquireSmallWorldsAtVertex(QTVertex*      givenVertex, rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world = rg_NULL);
    void inquireSmallWorldsAtEdge(  const QTEdge&  givenEdge,   rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world = rg_NULL);
    void inquireSmallWorldsAtFace(  const QTFace&  givenFace,   rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world = rg_NULL);
    void inquireSmallWorldsAtCell(  QTTetrahedron* givenCell,   rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world = rg_NULL);

    void inquireWholeSmallWorldsAtVertex(QTVertex* givenVertex, rg_dList<QTTetrahedron*>& smallWorldList, QTTetrahedron* world = rg_NULL);

    //  3. Query primitives for whole-world.
    //
    //  3.1. Q(v, Y): Query primitives for whole-world.
    void inquireVerticesOfVertex(QTVertex* givenVertex, rg_dList<QTVertex*>& vertexList           );
    void inquireEdgesOfVertex(   QTVertex* givenVertex, rg_dList<QTEdge>& edgeList               );
    void inquireFacesOfVertex(   QTVertex* givenVertex, rg_dList<QTFace>& faceList               );
    void inquireCellsOfVertex(   QTVertex* givenVertex, rg_dList<QTTetrahedron*>& tetrahedronList );
    
    //  3.2. Q(e, Y): Query primitives for whole-world.
    void inquireVerticesOfEdge(const QTEdge& givenEdge, rg_dList<QTVertex*>& vertexList           );
    void inquireEdgesOfEdge(   const QTEdge& givenEdge, rg_dList<QTEdge>& edgeList               );
    void inquireFacesOfEdge(   const QTEdge& givenEdge, rg_dList<QTFace>& faceList               );
    void inquireCellsOfEdge(   const QTEdge& givenEdge, rg_dList<QTTetrahedron*>& tetrahedronList );
    
    //  3.3. Q(f, Y): Query primitives for whole-world.
    void inquireVerticesOfFace(const QTFace& givenFace, rg_dList<QTVertex*>& vertexList           );
    void inquireEdgesOfFace(   const QTFace& givenFace, rg_dList<QTEdge>& edgeList               );
    void inquireFacesOfFace(   const QTFace& givenFace, rg_dList<QTFace>& faceList               );
    void inquireCellsOfFace(   const QTFace& givenFace, rg_dList<QTTetrahedron*>& tetrahedronList );
    
    //  3.4. Q(c, Y): Query primitives for whole-world.
    void inquireVerticesOfCell(QTTetrahedron* givenCell, rg_dList<QTVertex*>& vertexList           );
    void inquireEdgesOfCell(   QTTetrahedron* givenCell, rg_dList<QTEdge>& edgeList               );
    void inquireFacesOfCell(   QTTetrahedron* givenCell, rg_dList<QTFace>& faceList               );
    void inquireCellsOfCell(   QTTetrahedron* givenCell, rg_dList<QTTetrahedron*>& tetrahedronList );

    //
    //  END of Quasi-operator
    //  
    ///////////////////////////////////////////////////////////////////////////


    Sphere computeMaxTangentSphereOnQTTetraheron(QTTetrahedron const* ptrTetrahedron);
    Sphere computeMinSphereTangentToTwoBalls(const Sphere& ball1, const Sphere& ball2);
    rg_INT computeMinSphereTangentToThreeBalls(const Sphere& ball1, const Sphere& ball2, const Sphere& ball3, Sphere* tangentSphere);
        rg_INT makeCircumcircleOnCenterPlane(const Sphere& s1,
											 const Sphere& s2,
											 const Sphere& s3,
											 Sphere& result1,
											 Sphere& result2);
	        rg_INT makeCircumcircle(const rg_Circle2D& circle1,
							        const rg_Circle2D& circle2,
							        const rg_Circle2D& circle3,
							        rg_Circle2D&	result1, 
							        rg_Circle2D&	result2);
	        void shrinkCircle(const		rg_Circle2D& c1,
					          const		rg_Circle2D& c2,
					          const		rg_Circle2D& c3,
					          rg_Circle2D&	c1Tilde,
					          rg_Circle2D&	c2Tilde,
					          rg_Point2D&		c3Tilde);
	        rg_Point2D makeTangentLinesInWPlane(const rg_Circle2D&		c1, 
										        const rg_Circle2D&		c2,
										        const rg_Circle2D&		c3,
										        rg_ImplicitEquation&		result1,
										        rg_ImplicitEquation&		result2);

	        void makeTangentLinesOf2Circles(const rg_Circle2D&	circle1,
									        const rg_Circle2D&	circle2,
									        rg_ImplicitEquation&	result1,
									        rg_ImplicitEquation&	result2);

	        rg_Circle2D transformW2Z(const rg_ImplicitEquation& line,
							         const rg_Point2D& z3);
	        rg_Circle2D transformZ2W(const rg_Circle2D&	cTilde,
							         const rg_Point2D&	c3Tilde);

    void reportQuasiTriangulation( ofstream& fout );

    void analyzeAnomaly(ofstream& fout);


    void readQuasiTriangulationFile( const string& qtfFilename );
        int  changeLoadStateInQT(int& loadStateInQT, char* oneLineOfQTF);
        void passWhiteSpace(char* source, int& pos);
        void getNextSubstring(char* source, int& pos, char* destination);
        void loadQTVertex(QTVertex** qtVertex, char* oneLineOfQTF);
        void loadQTTetrahedron(QTVertex** qtVertex, QTTetrahedron** qtTetrahedron, char* oneLineOfQTF);
        void loadQTGate(QTVertex** qtVertex, QTTetrahedron** qtTetrahedron, char* oneLineOfQTF);
    void writeQuasiTriangulationFile( ofstream& fout, const string& source );

};

} // namespace GeometryTier

} // namespace V


#endif
