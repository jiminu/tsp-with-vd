#include "VertexGenerator2D.h"
using namespace V::GeometryTier;



VertexGenerator2D::VertexGenerator2D()
    : Generator2D(VERTEX_G)
    , m_childDiskGenerator(NULL)
    , m_vertexType(UNKNOWN_VTYPE)
    , m_prevEdgeGenerator(NULL)
    , m_nextEdgeGenerator(NULL)
    , m_isThisFromContainer(false)
{
}



VertexGenerator2D::VertexGenerator2D( const rg_Point2D & point, const int & userID )
    : Generator2D(VERTEX_G, userID)
    , m_point(point)
    , m_childDiskGenerator(NULL)
    , m_vertexType(UNKNOWN_VTYPE)
    , m_prevEdgeGenerator(NULL)
    , m_nextEdgeGenerator(NULL)
    , m_isThisFromContainer(false)
{
}



VertexGenerator2D::VertexGenerator2D( const VertexGenerator2D & vertexGenerator )
    : Generator2D(vertexGenerator)
    , m_point(vertexGenerator.m_point)
    , m_childDisk(vertexGenerator.m_childDisk)
    , m_childDiskGenerator(vertexGenerator.m_childDiskGenerator)
    , m_vertexType(vertexGenerator.m_vertexType)
    , m_prevEdgeGenerator(vertexGenerator.m_prevEdgeGenerator)
    , m_nextEdgeGenerator(vertexGenerator.m_nextEdgeGenerator)
    , m_isThisFromContainer(vertexGenerator.m_isThisFromContainer)
{
}



VertexGenerator2D::~VertexGenerator2D()
{
}



VertexGenerator2D& VertexGenerator2D::operator=( const VertexGenerator2D& generator )
{
    if (this == &generator)
        return *this;

    m_type                  = generator.m_type;
    m_ID                    = generator.m_ID;
    m_disk                  = generator.m_disk;
    m_outerFace             = generator.m_outerFace;
    m_innerFace             = generator.m_innerFace;
    m_innerVD               = generator.m_innerVD;
    m_userData              = generator.m_userData;
    m_innerGens             = generator.m_innerGens;

    m_point                 = generator.m_point;
    m_childDisk             = generator.m_childDisk;
    m_childDiskGenerator    = generator.m_childDiskGenerator;
    m_vertexType            = generator.m_vertexType;
    m_prevEdgeGenerator     = generator.m_prevEdgeGenerator;
    m_nextEdgeGenerator     = generator.m_nextEdgeGenerator;

    m_isThisFromContainer   = generator.m_isThisFromContainer;

    return *this;
}



bool VertexGenerator2D::operator==( const VertexGenerator2D & generator )
{
    if (this == &generator)
        return true;

    if (   m_type               == generator.m_type
        && m_ID                 == generator.m_ID
        && m_outerFace          == generator.m_outerFace
        && m_innerFace          == generator.m_innerFace
        && m_innerVD            == generator.m_innerVD
        && m_userData           == generator.m_userData
        && m_point              == generator.m_point
        && m_childDisk          == generator.m_childDisk
        && m_childDiskGenerator == generator.m_childDiskGenerator
        && m_vertexType         == generator.m_vertexType
        && m_prevEdgeGenerator  == generator.m_prevEdgeGenerator
        && m_nextEdgeGenerator  == generator.m_nextEdgeGenerator
        && m_isThisFromContainer== generator.m_isThisFromContainer)
    {
        return true;
    }
    else {
        return false;
    }
}


