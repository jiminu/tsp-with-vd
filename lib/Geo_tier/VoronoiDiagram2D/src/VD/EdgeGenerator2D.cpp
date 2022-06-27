#include "EdgeGenerator2D.h"
using namespace V::GeometryTier;



EdgeGenerator2D::EdgeGenerator2D()
    : Generator2D(EDGE_G)
    , m_firstSonDiskGenerator(NULL)
    , m_startVertexGenerator(NULL)
    , m_endVertexGenerator(NULL)
    , m_isThisFromContainer(false)
{
}



EdgeGenerator2D::EdgeGenerator2D( const rg_Point2D & startPoint, const rg_Point2D & endPoint, const rg_INT & userID )
    : Generator2D(EDGE_G, userID)
    , m_startPoint(startPoint)
    , m_endPoint(endPoint)
    , m_firstSonDiskGenerator(NULL)
    , m_startVertexGenerator(NULL)
    , m_endVertexGenerator(NULL)
    , m_isThisFromContainer(false)
{
}



EdgeGenerator2D::EdgeGenerator2D( const EdgeGenerator2D & edgeGenerator )
    : Generator2D(              edgeGenerator)
    , m_startPoint(             edgeGenerator.m_startPoint)
    , m_endPoint(               edgeGenerator.m_endPoint)
    , m_firstSonDiskGenerator(  edgeGenerator.m_firstSonDiskGenerator)
    , m_startVertexGenerator(   edgeGenerator.m_startVertexGenerator)
    , m_endVertexGenerator(     edgeGenerator.m_endVertexGenerator)
    , m_isThisFromContainer(    edgeGenerator.m_isThisFromContainer )
{
    m_childrenDisks             = edgeGenerator.m_childrenDisks;
    m_childrenDiskGenerators    = edgeGenerator.m_childrenDiskGenerators;
}



EdgeGenerator2D::~EdgeGenerator2D()
{
}



EdgeGenerator2D& EdgeGenerator2D::operator=( const EdgeGenerator2D& generator )
{
    if (this == &generator)
        return *this;

    m_type                      = generator.m_type;
    m_ID                        = generator.m_ID;
    m_disk                      = generator.m_disk;
    m_outerFace                 = generator.m_outerFace;
    m_innerFace                 = generator.m_innerFace;
    m_innerVD                   = generator.m_innerVD;
    m_userData                  = generator.m_userData;
    m_innerGens                 = generator.m_innerGens;

    m_startPoint                = generator.m_startPoint;
    m_endPoint                  = generator.m_endPoint;
    m_childrenDisks             = generator.m_childrenDisks;
    m_childrenDiskGenerators    = generator.m_childrenDiskGenerators;
    m_firstSonDiskGenerator     = generator.m_firstSonDiskGenerator;
    m_startVertexGenerator      = generator.m_startVertexGenerator;
    m_endVertexGenerator        = generator.m_endVertexGenerator;

    m_isThisFromContainer       = generator.m_isThisFromContainer;

    return *this;
}



bool EdgeGenerator2D::operator==( const EdgeGenerator2D & generator )
{
    if (this == &generator)
        return true;

    if (   m_type                   == generator.m_type
        && m_ID                     == generator.m_ID
        && m_outerFace              == generator.m_outerFace
        && m_innerFace              == generator.m_innerFace
        && m_innerVD                == generator.m_innerVD
        && m_userData               == generator.m_userData
        && m_startPoint             == generator.m_startPoint
        && m_endPoint               == generator.m_endPoint
        && m_firstSonDiskGenerator  == generator.m_firstSonDiskGenerator
        && m_startVertexGenerator   == generator.m_startVertexGenerator
        && m_endVertexGenerator     == generator.m_endVertexGenerator
        && m_isThisFromContainer    == generator.m_isThisFromContainer)
    {
        return true;
    }
    else {
        return false;
    }
}
