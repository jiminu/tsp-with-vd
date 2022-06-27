#ifndef _PATHINVESTIGATOR_H
#define _PATHINVESTIGATOR_H

#include "rg_Const.h"

#include "rg_Point3D.h"
#include "rg_dList.h"

#include "VDEdge.h"
#include "VDVertex.h"

namespace V {

namespace GeometryTier {


class PathInvestigator
{
private:
    VDVertex* m_vertex;
    VDEdge*   m_paths[3];
    rg_FLAG   m_orientation[3];
    rg_INT    m_currPathID;

public:
    PathInvestigator();
    PathInvestigator(const PathInvestigator& anInvestigator);
    ~PathInvestigator();

    VDVertex* getPropagatingVertex() const;
    VDEdge*  getPath(const rg_INT& i) const;
    VDEdge** getPaths();
    rg_FLAG  getOrientation(const rg_INT& i) const;
    VDEdge*  getCurrPath() const;
    rg_INT   getCurrPathID() const;
    rg_FLAG   getOrientationOfCurrPath() const;

    VDVertex* getNextPropagatingVertex() const;

    rg_FLAG  isValidInvestigator() const;

    void setPropagatingVertex(VDVertex* pVertex);
    void setPath(const rg_INT& i, VDEdge* aPath, const rg_FLAG& orientation);
//    void setPaths(VDEdge** paths);
//    void setCurrPathID(const rg_INT& id);
//    void setOrientation(const rg_FLAG& orientation);
    void setCurrPath(const rg_INT& id);
//    void setCurrPath(const rg_INT& id, VDEdge* aPath, const rg_FLAG& orientation);

    rg_FLAG makeNextInvestigator(const rg_Point3D& targetPt, rg_dList<VDVertex*>& verticesOnTrajectory, PathInvestigator& nextInvestigator);
    void    killCurrPath();

    PathInvestigator& operator =(const PathInvestigator& anInvestigator);
};

} // namespace GeometryTier

} // namespace V


#endif

