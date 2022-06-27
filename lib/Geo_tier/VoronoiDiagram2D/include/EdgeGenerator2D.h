#ifndef EDGE_GENERATOR2D_H
#define EDGE_GENERATOR2D_H

#include "rg_Circle2D.h"
#include "Generator2D.h"
#include "rg_Line2D.h"

namespace V {
namespace GeometryTier {

class VertexGenerator2D;
class DiskGenerator2D;

class EdgeGenerator2D : public Generator2D
{
private:
    rg_Point2D              m_startPoint;
    rg_Point2D              m_endPoint;
    list<rg_Circle2D>       m_childrenDisks;
    list<Generator2D*>      m_childrenDiskGenerators;
    Generator2D*            m_firstSonDiskGenerator;

    VertexGenerator2D*      m_startVertexGenerator;
    VertexGenerator2D*      m_endVertexGenerator;

    bool                    m_isThisFromContainer;

public:
    EdgeGenerator2D();
    EdgeGenerator2D(const rg_Point2D& startPoint, const rg_Point2D& endPoint, const rg_INT& userID);
    EdgeGenerator2D(const EdgeGenerator2D& edgeGenerator);
    ~EdgeGenerator2D();
    
    rg_Point2D              get_start_vertex_point() const;
    rg_Point2D              get_end_vertex_point() const;
    void                    get_geometry(rg_Line2D& lineSeg) const;
    rg_Line2D               get_geometry() const;
    list<rg_Circle2D>&      get_children_disks();
    list<Generator2D*>&     get_children_disk_generators();
    Generator2D*            get_first_son_disk_generator() const;
    VertexGenerator2D*      get_start_vertex_generator() const;
    VertexGenerator2D*      get_end_vertex_generator() const;
    

    void                    set_start_vertex_point(const rg_Point2D& startPt);
    void                    set_end_vertex_point(const rg_Point2D& endPt);
    void                    set_children_disks(const list<rg_Circle2D>& childrenDisks);
    void                    set_children_disk_generators(const list<Generator2D*>& childrenDiskGenerators);
    void                    set_first_son_disk_generator(const Generator2D* firstSonDiskGenerator);
    void                    set_start_vertex_generator(const VertexGenerator2D* startVtxGen);
    void                    set_end_vertex_generator(const VertexGenerator2D* endVtxGen);
    
    EdgeGenerator2D&        operator=(  const EdgeGenerator2D& generator);
    bool                    operator==( const EdgeGenerator2D& generator);

    bool                    isThisFromContainer() const;
    void                    isThisFromContainer( const bool& true_or_false );
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// inline functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline rg_Point2D EdgeGenerator2D::get_start_vertex_point() const { return m_startPoint; }
inline rg_Point2D EdgeGenerator2D::get_end_vertex_point() const { return m_endPoint; }
inline void EdgeGenerator2D::get_geometry(rg_Line2D& lineSeg) const { lineSeg.setLine2D(m_startPoint, m_endPoint); }
inline rg_Line2D EdgeGenerator2D::get_geometry() const { return rg_Line2D(m_startPoint, m_endPoint); }
inline VertexGenerator2D* EdgeGenerator2D::get_start_vertex_generator() const { return m_startVertexGenerator; }
inline VertexGenerator2D* EdgeGenerator2D::get_end_vertex_generator() const { return m_endVertexGenerator; }
inline list<Generator2D*>& EdgeGenerator2D::get_children_disk_generators() { return m_childrenDiskGenerators; }
inline Generator2D* EdgeGenerator2D::get_first_son_disk_generator() const { return m_firstSonDiskGenerator; }
inline list<rg_Circle2D>& EdgeGenerator2D::get_children_disks() { return m_childrenDisks; }

inline void EdgeGenerator2D::set_start_vertex_point(const rg_Point2D& startPt) { m_startPoint = startPt; }
inline void EdgeGenerator2D::set_end_vertex_point(const rg_Point2D& endPt) { m_endPoint = endPt; }
inline void EdgeGenerator2D::set_children_disks(const list<rg_Circle2D>& childrenDisks) { m_childrenDisks = childrenDisks; }
inline void EdgeGenerator2D::set_children_disk_generators(const list<Generator2D*>& childrenDiskGenerators) { m_childrenDiskGenerators = childrenDiskGenerators; }
inline void EdgeGenerator2D::set_first_son_disk_generator(const Generator2D* firstSonDiskGenerator) { m_firstSonDiskGenerator = const_cast<Generator2D*>(firstSonDiskGenerator); }
inline void EdgeGenerator2D::set_start_vertex_generator(const VertexGenerator2D* startVtxGen) { m_startVertexGenerator = const_cast<VertexGenerator2D*>(startVtxGen); }
inline void EdgeGenerator2D::set_end_vertex_generator(const VertexGenerator2D* endVtxGen) { m_endVertexGenerator = const_cast<VertexGenerator2D*>(endVtxGen); }

inline bool             EdgeGenerator2D::isThisFromContainer() const { return m_isThisFromContainer; }
inline void             EdgeGenerator2D::isThisFromContainer( const bool& true_or_false ) { m_isThisFromContainer = true_or_false; }

} // GeometryTier
} // V

#endif

