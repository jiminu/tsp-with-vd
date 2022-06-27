#include "FaceQT2D.h"
#include "VertexQT2D.h"

using namespace BULL2D::GeometryTier;


FaceQT2D::FaceQT2D()
{
}



FaceQT2D::FaceQT2D(const rg_INT& ID, const rg_Circle2D& tangentCircle)
: FaceSCDS(ID),
  m_emptyTangentCircle(tangentCircle)
{
}



FaceQT2D::FaceQT2D(const FaceQT2D& face)
: FaceSCDS(face),
  m_emptyTangentCircle(face.m_emptyTangentCircle)
{
}



FaceQT2D::~FaceQT2D()
{
}




rg_Circle2D FaceQT2D::getEmptyTangentCircle() const
{
    return m_emptyTangentCircle;
}



rg_BOOL FaceQT2D::isVirtual() const
{
    for ( rg_INT i=0; i<3; i++ ) {
        if ( ((VertexQT2D*)m_vertex[i])->isVirtual() ) {
            return rg_TRUE;
        }
    }

    return rg_FALSE;
}



void FaceQT2D::setEmptyTangentCircle(const rg_Circle2D& tangentCircle)
{
    m_emptyTangentCircle = tangentCircle;
}




FaceQT2D& FaceQT2D::operator =(const FaceQT2D& face)
{
    if ( this == &face ) {
        return *this;
    }

    FaceSCDS::operator =(face);

    m_emptyTangentCircle = face.m_emptyTangentCircle;

    return *this;
}



///////////////////////////////////////////////////////////////////////////

rg_BOOL FaceQT2D::isAnomaly() const
{
    if ( m_neighbor[0]==m_neighbor[1] || m_neighbor[0]==m_neighbor[2] || m_neighbor[1]==m_neighbor[2] ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}
