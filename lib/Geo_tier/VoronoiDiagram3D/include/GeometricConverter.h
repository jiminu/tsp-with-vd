#ifndef _GEOMETRICCONVERTER_H
#define _GEOMETRICCONVERTER_H


#include "rg_SphereSetVoronoiDiagram.h"
#include "rg_QuasiTriangulation.h"
#include "BetaUniverse.h"


namespace V {

namespace GeometryTier {


class GeometricConverter
{
private:
    GeometricConverter();

public:
    ~GeometricConverter();

    //static convertVoronoiDiagramIntoQuasiTriangulation(const rg_SphereSetVoronoiDiagram& VD, QuasiTriangulation& QT);
    //static convertQuasiTriangulationIntoVoronoiDiagram(const QuasiTriangulation& QT, rg_SphereSetVoronoiDiagram& VD, );
    static void convertQuasiTriangulationIntoVoronoiDiagram(const BetaUniverse& QTIneIWDS, rg_SphereSetVoronoiDiagram& VD);

private:
    //  Quasi-triangulation in eIWDS --> Voronoi diagram of spheres
    static void makeVVerticesCorrespondingToQCells( 
                                             const BetaUniverse&         QTIneIWDS,
                                             map<BetaCell*, VDVertex*>&  mapCellInBUToVertexInVD,
                                             rg_SphereSetVoronoiDiagram& VD );
    static void makeVEdgesCorrespondingToQFaces( 
                                             const BetaUniverse&         QTIneIWDS,
                                             map<BetaFace*, VDEdge*>&    mapFaceInBUToEdgeInVD,
                                             rg_SphereSetVoronoiDiagram& VD );
    static void makeVFacesCorrespondingToQEdges( 
                                             const BetaUniverse&         QTIneIWDS,
                                             map<BetaEdge*, VDFace*>&    mapEdgeInBUToFaceInVD,    
                                             rg_SphereSetVoronoiDiagram& VD );
    static void makeVCellsCorrespondingToQVertices( 
                                             const BetaUniverse&         QTIneIWDS,
                                             map<BetaVertex*, VDCell*>&  mapVertexInBUToCellInVD,
                                             rg_SphereSetVoronoiDiagram& VD );

    static void setTopologyOfVVerticesViaQCells( 
                                            const map<BetaCell*, VDVertex*>&  mapCellInBUToVertexInVD,
                                            const map<BetaFace*, VDEdge*>&    mapFaceInBUToEdgeInVD);
    static void setTopologyOfVEdgesViaQFaces( 
                                            const map<BetaFace*, VDEdge*>&    mapFaceInBUToEdgeInVD,
                                            const map<BetaCell*, VDVertex*>&  mapCellInBUToVertexInVD,
                                            const map<BetaEdge*, VDFace*>&    mapEdgeInBUToFaceInVD );
    static void setTopologyOfVFacesViaQEdges( 
                                            const map<BetaEdge*, VDFace*>&    mapEdgeInBUToFaceInVD,
                                            const map<BetaFace*, VDEdge*>&    mapFaceInBUToEdgeInVD,
                                            const map<BetaVertex*, VDCell*>&  mapVertexInBUToCellInVD );
    static void setTopologyOfVCellsViaQVertices( 
                                            const map<BetaVertex*, VDCell*>&  mapVertexInBUToCellInVD,
                                            const map<BetaEdge*, VDFace*>&    mapEdgeInBUToFaceInVD,
                                            rg_SphereSetVoronoiDiagram& VD );

};

} // namespace GeometryTier

} // namespace V


#endif

