#ifndef BULL2D_GEOMETRY_DELAUNAYTRIANGULATION2D_H
#define BULL2D_GEOMETRY_DELAUNAYTRIANGULATION2D_H

#include "rg_Const.h"
#include "rg_dList.h"

#include "Disc.h"
#include "VertexQT2D.h"
#include "FaceQT2D.h"

#include "VoronoiDiagram2DP.h"

#include <map>
using namespace std;

namespace BULL2D {
namespace GeometryTier {





class DelaunayTriangulation2D_GEO
{
private:
    //  quasi-simplicial complex data structure
    rg_dList<VertexQT2D> m_vertices;
    rg_dList<FaceQT2D>   m_faces;

    rg_dList<Disc>       m_discs;

public:
    DelaunayTriangulation2D_GEO();
    ~DelaunayTriangulation2D_GEO();

    inline rg_dList<VertexQT2D>& getVertices() {return m_vertices;};
    inline rg_dList<FaceQT2D>&   getFaces()    {return m_faces;};
    inline rg_dList<Disc>&       getBalls()    {return m_discs;};

	inline rg_INT getNumVertices() const       {return m_vertices.getSize();};
	inline rg_INT getNumFaces()    const       {return m_faces.getSize();};
	inline rg_INT getNumDiscs()    const       {return m_discs.getSize();};

    void construct(VoronoiDiagram2DP& pointVD);


    void getFacesIncidentToVertex(VertexQT2D* vertex, rg_dList<FaceQT2D*>& faceList) const;
    void getBoundingVerticesOfFace(FaceQT2D* face, rg_dList<VertexQT2D*>&  vertexList) const;
    void getNeighboringFacesOfFace(FaceQT2D* face, rg_dList<FaceQT2D*>&    faceList) const;
    void getBoundingVerticesOfFace(FaceQT2D* face, VertexQT2D** vertexArray) const;
    void getNeighboringFacesOfFace(FaceQT2D* face, FaceQT2D**   faceArray) const;

private:


	// 2014-04-15 JKKIM ADDED ///////////////////////////////////////////////////////////////////
	void createAndInitializeVerticesInQT2D(VoronoiDiagram2DP& pointSetVD, map<const VFace2DP*, VertexQT2D*>& mappingTableFromFaceInVD2VertexInQT);
	void createAndinitializeFacesInQT2D(   VoronoiDiagram2DP& pointSetVD, map<const VVertex2DP*, FaceQT2D*>& mappingTableFromVertexInVD2FaceInQT);
	void         setTopologyBetweenVerticesAndTetrahedraInQT(VoronoiDiagram2DP& pointSetVD,
															const map<const VFace2DP*, VertexQT2D*>& mappingTableFromFaceInVD2VertexInQT,
															const map<const VVertex2DP*, FaceQT2D*>&   mappingTableFromVertexInVD2FaceInQT);
	/////////////////////////////////////////////////////////////////////////////////////////////
};

} // namespace GeometryTier
} // namespace BULL2D

#endif


