#ifndef _QUASITRIANGULATION2D_H
#define _QUASITRIANGULATION2D_H

/******************************************************************************
 *    Related functions of CircleSetVoronoiDiagram and CircleVoronoiDiagram
 *    are commented out by CYSONG at July16, 2020.
 ******************************************************************************/

#include "rg_Const.h"
#include "rg_dList.h"
#include "Disc.h"
#include "VertexQT2D.h"
#include "FaceQT2D.h"
// #include "CircleVoronoiDiagram.h"
// #include "CircleSetVoronoiDiagram.h"
#include "VoronoiDiagram2DC.h"
#include <map>
using namespace std;


namespace V {
namespace GeometryTier {


class QuasiTriangulation2D
{
private:
    //  quasi-simplicial complex data structure
    rg_dList<VertexQT2D> m_vertices;
    rg_dList<FaceQT2D>   m_faces;

    rg_dList<Disc>       m_discs;

public:
    QuasiTriangulation2D();
    ~QuasiTriangulation2D();

    inline rg_dList<VertexQT2D>& getVertices() {return m_vertices;};
    inline rg_dList<FaceQT2D>&   getFaces()    {return m_faces;};
    inline rg_dList<Disc>&       getBalls()    {return m_discs;};

	inline rg_INT getNumVertices() const       {return m_vertices.getSize();};
	inline rg_INT getNumFaces()    const       {return m_faces.getSize();};
	inline rg_INT getNumDiscs()    const       {return m_discs.getSize();};

//     void construct(CircleVoronoiDiagram& circleVD);
//     void construct(CircleSetVoronoiDiagram& circleVD);
	void construct(VoronoiDiagram2DC& circleVD);


    void getFacesIncidentToVertex(VertexQT2D* vertex, rg_dList<FaceQT2D*>& faceList) const;
    void getBoundingVerticesOfFace(FaceQT2D* face, rg_dList<VertexQT2D*>&  vertexList) const;
    void getNeighboringFacesOfFace(FaceQT2D* face, rg_dList<FaceQT2D*>&    faceList) const;
    void getBoundingVerticesOfFace(FaceQT2D* face, VertexQT2D** vertexArray) const;
    void getNeighboringFacesOfFace(FaceQT2D* face, FaceQT2D**   faceArray) const;

private:
    /*
    VertexQT2D** createAndInitializeVerticesInQT2D(CircleVoronoiDiagram& circleSetVD);
    FaceQT2D**   createAndinitializeFacesInQT2D(   CircleVoronoiDiagram& circleSetVD);

    void setTopologyBetweenVerticesAndTetrahedraInQT(CircleVoronoiDiagram& circleSetVD,
                                                     VertexQT2D** mappingTableFromFaceInVD2VertexInQT,
                                                     FaceQT2D**   mappingTableFromVertexInVD2FaceInQT);


    VertexQT2D** createAndInitializeVerticesInQT2D(CircleSetVoronoiDiagram& circleSetVD);
    FaceQT2D**   createAndinitializeFacesInQT2D(   CircleSetVoronoiDiagram& circleSetVD);
    void         setTopologyBetweenVerticesAndTetrahedraInQT(CircleSetVoronoiDiagram& circleSetVD,
															VertexQT2D** mappingTableFromFaceInVD2VertexInQT,
															FaceQT2D**   mappingTableFromVertexInVD2FaceInQT);
    */

	// 2014-04-15 JKKIM ADDED ///////////////////////////////////////////////////////////////////
	void createAndInitializeVerticesInQT2D(VoronoiDiagram2DC& circleSetVD, map<const VFace2D*, VertexQT2D*>& mappingTableFromFaceInVD2VertexInQT);
	void createAndinitializeFacesInQT2D(   VoronoiDiagram2DC& circleSetVD, map<const VVertex2D*, FaceQT2D*>& mappingTableFromVertexInVD2FaceInQT);
	void         setTopologyBetweenVerticesAndTetrahedraInQT(VoronoiDiagram2DC& circleSetVD,
															const map<const VFace2D*, VertexQT2D*>& mappingTableFromFaceInVD2VertexInQT,
															const map<const VVertex2D*, FaceQT2D*>&   mappingTableFromVertexInVD2FaceInQT);
	/////////////////////////////////////////////////////////////////////////////////////////////
};


} // GeometryTier
} // V

#endif

