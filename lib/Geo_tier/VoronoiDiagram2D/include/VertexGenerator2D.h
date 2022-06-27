#ifndef VERTEX_GENERATOR2D_H
#define VERTEX_GENERATOR2D_H

#include "Generator2D.h"
#include "rg_Circle2D.h"

namespace V {
namespace GeometryTier {

class EdgeGenerator2D;
class DiskGenerator2D;

class VertexGenerator2D : public Generator2D
{
public:
    enum Vertex_Type { REFLEX_FROM_POLYGON_INSIDE, NON_REFLEX_FROM_POLYGON_INSIDE, UNKNOWN_VTYPE };

private:
    rg_Point2D          m_point;
    rg_Circle2D         m_childDisk;
    Generator2D*        m_childDiskGenerator;
    Vertex_Type         m_vertexType;

    EdgeGenerator2D*    m_prevEdgeGenerator;
    EdgeGenerator2D*    m_nextEdgeGenerator;

    bool                m_isThisFromContainer;


public:
    VertexGenerator2D();
    VertexGenerator2D(const rg_Point2D& point, const int& userID);
    VertexGenerator2D(const VertexGenerator2D& vertexGenerator);
    ~VertexGenerator2D();
    
    
    rg_Point2D          get_point() const;
    void                get_geometry(rg_Point2D& point) const;
    rg_Circle2D         get_child_disk() const;
    Generator2D*        get_child_disk_generator() const;
    rg_Circle2D*        get_child_disk_ptr();
    Generator2D*        get_first_son_disk_generator()  const;
    Vertex_Type         vertex_type() const;
    EdgeGenerator2D*    get_next_edge_generator() const;
    EdgeGenerator2D*    get_previous_edge_generator() const;
    

    void                set_point(const rg_Point2D& point);
    void                set_child_disk(const rg_Circle2D& childDisk);
    void                set_child_disk_generator(const Generator2D* childDiskGenerator);
    void                set_first_son_disk_generator(const Generator2D* firstSonDiskGenerator);
    void                set_vertex_type(const Vertex_Type& vertexType);
    void                set_next_edge_generator( const EdgeGenerator2D* nextEdgeGen );
    void                set_previous_edge_generator( const EdgeGenerator2D* prevEdgeGen );


    VertexGenerator2D&  operator=(  const VertexGenerator2D& generator);
    bool                operator==( const VertexGenerator2D& generator);
    //void             duplicate(const Generator2D& generator);

    bool                isThisFromContainer() const;
    void                isThisFromContainer( const bool& true_or_false );
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// inline functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline rg_Point2D       VertexGenerator2D::get_point() const { return m_point; }
inline void             VertexGenerator2D::get_geometry(rg_Point2D& point)  const { point = m_point; }
inline rg_Circle2D      VertexGenerator2D::get_child_disk() const { return m_childDisk; }
inline rg_Circle2D*     VertexGenerator2D::get_child_disk_ptr() { return &m_childDisk; }
inline Generator2D*     VertexGenerator2D::get_child_disk_generator() const { return m_childDiskGenerator; }
inline Generator2D*     VertexGenerator2D::get_first_son_disk_generator()  const { return m_childDiskGenerator;  }
inline VertexGenerator2D::Vertex_Type VertexGenerator2D::vertex_type() const { return m_vertexType; }
inline EdgeGenerator2D* VertexGenerator2D::get_next_edge_generator() const { return m_nextEdgeGenerator; }
inline EdgeGenerator2D* VertexGenerator2D::get_previous_edge_generator() const { return m_prevEdgeGenerator; }


inline void             VertexGenerator2D::set_point(const rg_Point2D& point) { m_point = point; }
inline void             VertexGenerator2D::set_child_disk(const rg_Circle2D& childDisk) { m_childDisk = childDisk; }
inline void             VertexGenerator2D::set_child_disk_generator(const Generator2D* childDiskGenerator) { m_childDiskGenerator = const_cast<Generator2D*>(childDiskGenerator); }
inline void             VertexGenerator2D::set_first_son_disk_generator(const Generator2D* firstSonDiskGenerator) {  m_childDiskGenerator = const_cast<Generator2D*>(firstSonDiskGenerator); }
inline void             VertexGenerator2D::set_vertex_type(const Vertex_Type& vertexType) { m_vertexType = vertexType; }
inline void             VertexGenerator2D::set_next_edge_generator( const EdgeGenerator2D* nextEdgeGen ) { m_nextEdgeGenerator = const_cast<EdgeGenerator2D*>(nextEdgeGen); }
inline void             VertexGenerator2D::set_previous_edge_generator( const EdgeGenerator2D* prevEdgeGen ) { m_prevEdgeGenerator = const_cast<EdgeGenerator2D*>(prevEdgeGen); }

inline bool             VertexGenerator2D::isThisFromContainer() const { return m_isThisFromContainer; }
inline void             VertexGenerator2D::isThisFromContainer( const bool& true_or_false ) { m_isThisFromContainer = true_or_false; }

} // GeometryTier
} // V

#endif
