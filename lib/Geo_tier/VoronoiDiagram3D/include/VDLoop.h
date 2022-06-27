#ifndef _VDLOOP_H
#define _VDLOOP_H

#include "ConstForVoronoiDiagram3D.h"

#include "TopologicalEntity.h"
#include "rg_dList.h"

namespace V {

namespace GeometryTier {


class VDFace;
class VDPartialEdge;
class VDEdge;

class VDLoop : public TopologicalEntity
{
private:
    VDFace*        m_face;
    VDPartialEdge* m_partEdge;

    rg_FLAG        m_isOuterLoop;

    rg_INT         m_numOfBoundingEdges;

public:
    //  constructor & deconstructor..
    VDLoop();
    VDLoop(const rg_INT& ID, VDFace* face);
    VDLoop(const rg_INT& ID, VDFace* face, const rg_FLAG& bIsOuterLoop );
    VDLoop(const rg_INT& ID, VDFace* face, VDPartialEdge* partEdge, const rg_FLAG& bIsOuterLoop);
    VDLoop(const TopologicalEntity& aTopoEntity, VDFace* face, VDPartialEdge* partEdge, const rg_FLAG& bIsOuterLoop);
    VDLoop(const VDLoop& aLoop);
    ~VDLoop();

    //  get functions.. 
    VDFace*        getFace() const;
    VDPartialEdge* getPartialEdge() const;
    rg_INT         getNumOfBoundingEdges() const;
    rg_FLAG        isOuterLoop() const;

    rg_FLAG        isClosedLoop() const;
    rg_BOOL        isSingleEdge() const;

    //  set functions..
    void setFace(VDFace* face);
    void setPartialEdge(VDPartialEdge* partEdge);
    void isOuterLoop(const rg_FLAG& bIsOuterLoop);

    void setLoop(VDFace* face, VDPartialEdge* partEdge, const rg_FLAG& bIsOuterLoop);
    void setLoop(const rg_INT& ID, VDFace* face, VDPartialEdge* partEdge, const rg_FLAG& bIsOuterLoop);

    rg_FLAG addPartialEdgeByEdgeTracing(VDPartialEdge *const partialEdgeToAdd);
    void addPartialEdgeToLoop(VDPartialEdge *const partEdgeToAdd);
    void addPartialEdgeOnInfinityToLoop(VDPartialEdge *const partEdgeToAdd);

    rg_FLAG defineOrderOfPartialEdgesInLoop( rg_dList<VDPartialEdge*>& multipleOuterLoop );
    rg_FLAG reorderOfPartialEdgesInLoop( rg_dList<VDPartialEdge*>& multipleOuterLoop );
    rg_FLAG constructLoopCycle( rg_dList<VDPartialEdge*>& multipleOuterLoop );

    //  operator overloading..
    VDLoop& operator =(const VDLoop& aLoop);

    //  topological operators..
    void collectBoundingPartialEdges( rg_dList<VDPartialEdge*>& boundingPartialEdges);
    void collectBoundingPartialEdgesInCCW( rg_dList<VDPartialEdge*>& boundingPartialEdges);
    void collectBoundingPartialEdgesInCW( rg_dList<VDPartialEdge*>& boundingPartialEdges);


    rg_BOOL searchBoundingEdges( rg_dList<VDEdge*>& boundingEdges ) const;
};

} // namespace GeometryTier

} // namespace V


#endif

