#ifndef GENERATOR2D_H
#define GENERATOR2D_H

#include "rg_Point2D.h"
#include "rg_Circle2D.h"



#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {


class VFace2D;

class Generator2D
{
private:
    int             m_ID;
    rg_Circle2D*    m_disk;
    VFace2D*        m_assignedVFace;


public:
    Generator2D();
    Generator2D( const int& ID );
    Generator2D( const int& ID, rg_Circle2D* const disk );
    Generator2D( const int& ID, rg_Circle2D* const disk, VFace2D* const face );
    Generator2D( const Generator2D& generator );
    ~Generator2D();

    inline int          getID() const { return m_ID; }
    inline rg_Circle2D* getDisk() const { return m_disk; }
    inline VFace2D*     getAssignedVFace() const { return m_assignedVFace; }

    void                setID( const int& ID ) { m_ID = ID; }
    void                setDisk( rg_Circle2D* const disk ) { m_disk = disk; }
    void                setAssignedVFace( VFace2D* const face ) { m_assignedVFace = face; }

    Generator2D&        operator=( const Generator2D& generator );
    bool                operator==( const Generator2D& generator ) const;

    void                getNeighborGenerators( list<Generator2D*>& neighborGeneratorsList );

};

} // namespace GeometryTier
} // namespace BULL2D

#endif


