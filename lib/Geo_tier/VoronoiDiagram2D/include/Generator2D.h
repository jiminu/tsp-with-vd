#ifndef GENERATOR2D_H
#define GENERATOR2D_H

#include "rg_Point2D.h"
#include "rg_Circle2D.h"

#include <list>
using namespace std;

namespace V {
namespace GeometryTier {

class VFace2D;
class VoronoiDiagram2DC;
class VoronoiDiagramCIC;

class Generator2D
{
public:
    enum Generator_Type { POINT_G, DISK_G, ELLIPSE_G, LINESEGMENT_G, SPLINE_G, EDGE_G, VERTEX_G, POLYGON_G, CHAIN_G, UNKNOWN_G };


protected:
    Generator_Type      m_type;
    int                 m_ID;
    rg_Circle2D         m_disk;
    VFace2D*            m_outerFace;
    VFace2D*            m_innerFace;
    VoronoiDiagramCIC*  m_innerVD;
    void*               m_userData;
    list<Generator2D*>  m_innerGens;


public:
    //BEGIN - previous constructor
    /* Many applications are using Generator2D as DiskGenerator2D 
    so we just let a Generator2D be a DiskGenerator2D and Diskgenerator2D.
    Previous constructors set m_type to DISK_G */
    Generator2D();
    Generator2D( const int& ID );
    Generator2D( const int& ID, const rg_Circle2D& disk );
    Generator2D( const int& ID, const rg_Circle2D& disk, VFace2D* const face );
    Generator2D( const int& ID, const rg_Circle2D& disk, void* const userData );
    Generator2D( const Generator2D& generator );
    //END - previous constructor

    //BEGIN - constructors for polygon VD
    Generator2D(const Generator_Type& type);
    Generator2D(const Generator_Type& type, const int& ID);
    Generator2D(const Generator_Type& type, const VFace2D* const Vface, void* userData=NULL, const int& ID=0);
    //BEGIN - constructors for polygon VD
    ~Generator2D();

    inline int                  getID() const               { return m_ID; }
    inline rg_Circle2D          getDisk() const             { return m_disk; }
    inline rg_Circle2D*         getDiskPtr()                { return &m_disk; } //CYSONG[Oct17, 19]: function added
	inline VFace2D*             getOuterFace() const        { return m_outerFace; }
    inline VFace2D*             getInnerFace() const        { return m_innerFace; }
    inline VoronoiDiagramCIC*   getInnerVD() const          { return m_innerVD; }
    inline void*                getUserData() const         { return m_userData; }
    inline void                 getInnerGens(list<Generator2D*>& innerGens) const { innerGens = m_innerGens; }
    inline Generator_Type       getType() const             { return m_type; }

    inline void                 setID( const int& ID )                      { m_ID = ID; }
    inline void                 setDisk( const rg_Circle2D& disk )          { m_disk = disk; }
    inline void                 setOuterFace( VFace2D* const face )         { m_outerFace = face; }
    inline void                 setInnerFace( VFace2D* const face )         { m_innerFace = face; }
    inline void                 setUserData( const void* const userData )   { m_userData = const_cast<void*>(userData); }
    inline void                 setInnerGens( const list<Generator2D*>& innerGens ) { m_innerGens = innerGens; }
    inline void                 setInnerVoronoiDiagram( VoronoiDiagramCIC* const innerVD ) { m_innerVD = innerVD; }
    inline void                 setType( const Generator_Type& type )       { m_type = type; }

    void                        addOneInnerGen( Generator2D* const innerGen ) { m_innerGens.push_back(innerGen); }
    bool                        constructVD();

    Generator2D&                operator=( const Generator2D& generator );
    bool                        operator==( const Generator2D& generator ) const;

    void                        getNeighborGeneratorsInThisVoronoiDiagram( list<Generator2D*>& neighborGeneratorsList, VoronoiDiagram2DC* const VD );
    void                        getNeighborGenerators(list<Generator2D*>& neighborGeneratorsList, const bool& isThisContainer) const;
    void                        removeOneInnerGen(Generator2D* const innerGen) { m_innerGens.remove(innerGen); };
    void                        clearInnerGen() { m_innerGens.clear(); };
};

} // GeometryTier
} // V

#endif


