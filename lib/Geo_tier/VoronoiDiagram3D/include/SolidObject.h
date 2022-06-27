#ifndef _SOLIDOBJECT_H
#define _SOLIDOBJECT_H

#include "rg_Const.h"
#include "rg_dList.h"
#include "rg_Point3D.h"

#include "rg_WingedEdgeDataStructure.h"

#include <fstream>
using namespace std;

namespace V {

namespace GeometryTier {


typedef WVertex< rg_Point3D, void*, void*> VertexOfSolid;
typedef WEdge<   rg_Point3D, void*, void*> EdgeOfSolid;
typedef WFace<   rg_Point3D, void*, void*> FaceOfSolid;

class SolidObject
{
private:
    rg_dList< VertexOfSolid > m_vertices;
    rg_dList< EdgeOfSolid >   m_edges;
    rg_dList< FaceOfSolid >   m_faces;

public:
    SolidObject();
    ~SolidObject();

    rg_dList< VertexOfSolid >* getVertices();
    rg_dList< EdgeOfSolid >*   getEdges();
    rg_dList< FaceOfSolid >*   getFaces();

    VertexOfSolid* splitEdge(EdgeOfSolid* edge, const rg_Point3D& splitPt);
    void makeEdge(FaceOfSolid* face, VertexOfSolid* startVertex, VertexOfSolid* endVertex);

    void makeUnitOctahedron();
    void makeSphericalTriangularMesh(const rg_INT& levelOfResolution, const rg_Point3D& center = rg_Point3D(0.0, 0.0, 0.0), const rg_REAL& radius = 1.0);
    void scaleSphericalTriangularMesh(const rg_REAL& radius);
    void translateSphericalTriangularMesh(const rg_Point3D& delta);

    void writeIndexedFacedSetOfVRML(char* firstIndention, char* indentionUnit, ofstream& fout);
    
};

} // namespace GeometryTier

} // namespace V


#endif
