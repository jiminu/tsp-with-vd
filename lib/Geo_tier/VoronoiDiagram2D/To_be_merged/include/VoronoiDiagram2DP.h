#ifndef VORONOIDIAGRAM2DP_H
#define VORONOIDIAGRAM2DP_H

#include "VEdge2DP.h"
#include "VFace2DP.h"
#include "VVertex2DP.h"
#include "rg_Point2D.h"
#include "Generator2DP.h"

#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {

class VoronoiDiagram2DP
{
private:
    list<VEdge2DP>       m_VEdges;
    list<VFace2DP>       m_VFaces;
    list<VVertex2DP>     m_VVertices;
    list<Generator2DP>   m_generators;
    list<Generator2DP>   m_phantomGenerators;

public:
    VoronoiDiagram2DP();
    VoronoiDiagram2DP( list<rg_Point2D>& pointset );
    VoronoiDiagram2DP( const VoronoiDiagram2DP& VD2DP );
    ~VoronoiDiagram2DP();

    void         getVoronoiEdges(    list<const VEdge2DP*>& VEdgesList ) const;
    void         getVoronoiFaces(    list<const VFace2DP*>& VFacesList ) const;
    void         getVoronoiVertices( list<const VVertex2DP*>& VVerticesList ) const;
    void         getGenerators(      list<const Generator2DP*>& generatorsList ) const;

    void         setGenerators(      list<rg_Point2D>& pointset );
    void         setGenerators(      list<Generator2DP>& generators );

    VoronoiDiagram2DP& operator=(    const VoronoiDiagram2DP& VD2DP );



    void         constructVoronoiDiagram( list<rg_Point2D>& pointset );
    void         constructVoronoiDiagram();
    VFace2DP*    findVFaceContainingQueryPoint( const rg_Point2D& pt );
    void         clear();

private:
    VVertex2DP*   createVVertex( const VVertex2DP& vertex );
    VEdge2DP*     createVEdge(   const VEdge2DP& edge );
    VFace2DP*     createVFace(   const VFace2DP& face );

    void          makePhantomGenerators();
    void          constructVoronoiDiagramForPhantomGenerators();
    void          updateVoronoiDiagramWithNewGenerator( Generator2DP* const newGenerator );
    Generator2DP* findClosestGeneratorToNewGenerator(   Generator2DP* const newGenerator );
    void          removeAllExtraneousVVerticesAndVEdges();

    void          findIntersectingVEdgesOfCurrVDAgainstNewVEdgeLoop( Generator2DP* const newGenerator,
                                                                     Generator2DP* const closestGenerator,
                                                                     list<VEdge2DP*>& intersectingVEdges );
    void          makeVEdgeLoopForNewGeneratorAndConnectToCurrVD(    Generator2DP* const newGenerator,
                                                                     list<VEdge2DP*>& intersectingVEdges,
                                                                     list<VVertex2DP*>& newVVertices,
                                                                     list<VEdge2DP*>& newVEdges );
    void          connectCurrVDToNewVEdgeLoop( list<VEdge2DP*>& newVEdges );
    void          computeCoordOfNewVVertices(  list<VVertex2DP*>& newVVertices );

    VVertex2DP*   findFirstRedVertex( Generator2DP* const newGenerator, Generator2DP* const closestGenerator );
    bool          hasCycleOccurred(   VVertex2DP* const candidateForRedVVertex );
    bool          hasOldRegionSplit(  VVertex2DP* const candidateForRedVVertex );

    rg_Point2D    computeCircumcenter( const rg_Point2D& pt1, const rg_Point2D& pt2, const rg_Point2D& pt3 );
};

} // namespace GeometryTier
} // namespace BULL2D


#endif


