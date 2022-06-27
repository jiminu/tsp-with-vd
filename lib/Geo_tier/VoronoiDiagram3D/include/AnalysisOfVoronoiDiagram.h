#ifndef _ANALYSISOFVORONOIDIAGRAM_H
#define _ANALYSISOFVORONOIDIAGRAM_H

#include "rg_SphereSetVoronoiDiagram.h"

#include <list>
using namespace std;


namespace V {

namespace GeometryTier {


class AnalysisOfVoronoiDiagram
{
private:
    rg_SphereSetVoronoiDiagram* m_VoronoiDiagram;


public:
    AnalysisOfVoronoiDiagram();
    AnalysisOfVoronoiDiagram(rg_SphereSetVoronoiDiagram* diagram);
    AnalysisOfVoronoiDiagram(const AnalysisOfVoronoiDiagram& analysis);
    ~AnalysisOfVoronoiDiagram();

    inline rg_SphereSetVoronoiDiagram* getVoronoiDiagram() { return m_VoronoiDiagram; }
    inline void                        setVoronoiDiagram(rg_SphereSetVoronoiDiagram* diagram) { m_VoronoiDiagram = diagram; }

    AnalysisOfVoronoiDiagram& operator =(const AnalysisOfVoronoiDiagram& analysis);

    void reportCorrectnessOfVDGeometry(ofstream& fout);


    void reportCorrectnessOfVertexGeometry(ofstream& fout);
    void reportCorrectnessOfEdgeGeometry(ofstream& fout);

};

} // namespace GeometryTier

} // namespace V


#endif

