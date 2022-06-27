#ifndef GENERATOR2DP_H
#define GENERATOR2DP_H

#include "rg_Point2D.h"



#include <list>
using namespace std;

namespace BULL2D {
namespace GeometryTier {


class VFace2DP;

class Generator2DP
{
private:
    int         m_ID;
    rg_Point2D* m_location;
    VFace2DP*   m_assignedVFace;

public:
    Generator2DP();
    Generator2DP( const int& ID );
    Generator2DP( const int& ID, rg_Point2D* const location );
    Generator2DP( const int& ID, rg_Point2D* const location, VFace2DP* const face );
    Generator2DP( const Generator2DP& generator );
    ~Generator2DP();

    inline int          getID() const { return m_ID; }
    inline rg_Point2D*  getLocation() const { return m_location; }
    inline VFace2DP*    getAssignedVFace() const { return m_assignedVFace; }

    void                setID( const int& ID ) { m_ID = ID; }
    void                setLocation( rg_Point2D* const location ) { m_location = location; }
    void                setAssignedVFace( VFace2DP* const face ) { m_assignedVFace = face; }

    Generator2DP&       operator=( const Generator2DP& generator );

    void                getNeighborGenerators( list<Generator2DP*>& neighborGeneratorsList );

};

} // namespace GeometryTier
} // namespace BULL2D

#endif


