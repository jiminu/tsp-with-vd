#ifndef TRIANGULATION3D_H
#define TRIANGULATION3D_H

#include "rg_dList.h"
#include "T3DVertex.h"
#include "T3DTetrahedron.h"
#include "T3DGate.h"
#include "IntegerBy96BIT.h"
#include "GateForT3D.h"

#include <fstream>
using namespace std;

namespace V {

namespace GeometryTier {


const rg_INT NO_INIT_CONDITION        = 0;
const rg_INT GIVEN_FACE_ON_CONVEXHULL = 1;

class Triangulation3D
{
private:
    rg_dList<T3DVertex>      m_vertices;
    rg_dList<T3DTetrahedron> m_tetrahedra;

    rg_INT m_decimalPoint;

//    ofstream FOUT;
public:
    Triangulation3D();
    ~Triangulation3D();

    void reset();

    rg_INT getNumOfVertices() const;
    rg_INT getNumOfTetrahedra() const;

    rg_dList<T3DVertex>*      getVertices();
    rg_dList<T3DTetrahedron>* getTetrahedra();

    T3DVertex*      addVertex(const T3DVertex& vertex);
    T3DTetrahedron* addTetrahedron(const T3DTetrahedron& tetrahedron);
    rg_FLAG         locateFaceByVertexID(T3DFace& faceToLocate, const rg_INT& vertexID1, const rg_INT& vertexID2, const rg_INT& vertexID3);

    void tetrahedralizePointsOnSphere(const rg_FLAG& initCondition, const rg_INT& numPts, rg_Point3D* points);
        void loadPoints(const rg_INT& numPts, rg_Point3D* points);
		void initializeTetrahedralizationOfPointsOnSphere(const rg_FLAG& initCondition, 
                                                          T3DVertex*& convergingVertex,
                                                          rg_dList< rg_dList<GateForT3D>* >& listOfMakingLine);
            rg_dList<GateForT3D>* makeInitialMakingLine(T3DVertex*& convergingVertex, T3DVertex** vertexForInitMakingLine);
		T3DVertex* findVertexToDefineAdjacentFaceOnConvexHull(T3DVertex* vertex1, T3DVertex* vertex2, T3DVertex* mateVertex);
            rg_FLAG do3VerticesDefineConvexHullFace(T3DVertex* vertex1, T3DVertex* vertex2, T3DVertex* vertex3);
                rg_FLAG isLastPointInPlusHalfSpaceByFirstThreePoints(rg_Point3D* point, const rg_Point3D& lastPoint);
                IntegerBy96BIT doConvexHullTestByExactArithmeticBasedOnSymbolicPerburbation(
                                   rg_INT* intPt1, rg_INT* intPt2, rg_INT* intPt3, rg_INT* intPt4);
        void updateMakingLineByNewVertex(T3DVertex* convergingVertex, 
                                            T3DVertex* vertexForAdjacentFace, 
                                            rg_dNode<GateForT3D>*& currNode,
                                            rg_dList<GateForT3D>* makingLine);
        void updateMakingLineByVertexOnLine(T3DVertex* convergingVertex, 
                                            T3DVertex* vertexForAdjacentFace, 
                                            rg_dNode<GateForT3D>*& currNode,
                                            rg_dList<GateForT3D>* makingLine,
                                            rg_dList< rg_dList<GateForT3D>* >& listOfMakingLine);
            
            void reportMakingLine(rg_dList<GateForT3D>* makingLine);



    void triangulatePointsGivenFaceOnConvexHull(const rg_INT& numPts, rg_Point3D* points);
		void initializeTriangulationGivenFaceOnConvexHull(rg_dList<T3DGate>& gateQueue);
		void findVertexToDefineConvexHullFaceWithEdgeOfCurrFace(T3DGate& currGate);

        void selectVertexToMakeConvexHullFace(T3DGate& currGate, rg_dList<T3DGate>& gateQueue);
            T3DTetrahedron* makeTetrahedron(T3DGate& currGate, const rg_INT& edgePos);
            T3DGate*        makeGate(rg_dList<T3DGate>& gateQueue, T3DTetrahedron* tetrahedron, 
                                     const rg_INT& mateVertexPos, const rg_FLAG& priority);
            T3DGate* isThereSameGate(rg_dList<T3DGate>& gateQueue, const T3DGate& newGate);
            T3DGate* isThereSameGate(rg_dList<T3DGate>& gateQueue, T3DVertex** vertexOfGate);
            rg_INT   exploreConnectionBetweenGatesAndVertices(rg_dList<T3DGate>& gateQueue, T3DGate& currGate);

    void loadRandomData(const rg_INT& numVertices);
    void loadTestData();

    void testDiscrimentFunction();
        rg_REAL computeDetermineOriginalFloatingPointArithmetic(const rg_Point3D& pt1, const rg_Point3D& pt2, 
                                                            const rg_Point3D& pt3, const rg_Point3D& pt4);
        rg_REAL computeDetermineFullFloatingPointArithmetic(const rg_Point3D& pt1, const rg_Point3D& pt2, 
                                                            const rg_Point3D& pt3, const rg_Point3D& pt4);
        rg_REAL computeDetermineFiniteFloatingPointArithmetic(rg_REAL& errorUpperBound, 
                                                              const rg_Point3D& pt1, const rg_Point3D& pt2, 
                                                              const rg_Point3D& pt3, const rg_Point3D& pt4);
};

} // namespace GeometryTier

} // namespace V


#endif
