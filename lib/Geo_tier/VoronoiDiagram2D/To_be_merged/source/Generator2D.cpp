#include "Generator2D.h"
#include "VFace2D.h"
using namespace BULL2D::GeometryTier;



Generator2D::Generator2D()
{
    m_ID            = -1;
    m_disk          = NULL;
    m_assignedVFace = NULL;
}



Generator2D::Generator2D( const int& ID )
{
    m_ID            = ID;
    m_disk          = NULL;
    m_assignedVFace = NULL;
}



Generator2D::Generator2D( const int& ID, rg_Circle2D* const disk )
{
    m_ID            = ID;
    m_disk          = disk;
    m_assignedVFace = NULL;
}



Generator2D::Generator2D( const int& ID, rg_Circle2D* const disk, VFace2D* const face )
{
    m_ID            = ID;
    m_disk          = disk;
    m_assignedVFace = face;
}



Generator2D::Generator2D( const Generator2D& generator )
{
    m_ID            = generator.m_ID;
    m_disk          = generator.m_disk;
    m_assignedVFace = generator.m_assignedVFace;
}



Generator2D::~Generator2D()
{
}



Generator2D& Generator2D::operator=( const Generator2D& generator )
{
    if( this == &generator ) {
        return *this;
    }

    m_ID            = generator.m_ID;
    m_disk          = generator.m_disk;
    m_assignedVFace = generator.m_assignedVFace;

    return *this;
}




bool Generator2D::operator==( const Generator2D& generator ) const
{
    if( this == &generator ) {
        return true;
    }

    if( m_ID == generator.m_ID &&
        m_disk == generator.m_disk &&
        m_assignedVFace == generator.m_assignedVFace ) 
    {
        return true;
    }
    else {
        return false;
    }
}



void Generator2D::getNeighborGenerators( list<Generator2D*>& neighborGeneratorsList )
{
    list<VFace2D*> adjacentVFaces;
    m_assignedVFace->getAdjacentVFaces( adjacentVFaces );

    for( list<VFace2D*>::iterator it_face = adjacentVFaces.begin() ; it_face != adjacentVFaces.end() ; it_face++ ) {
        VFace2D* currVFace = *it_face;

        if( currVFace->getGenerator() != NULL ) {
            neighborGeneratorsList.push_back( currVFace->getGenerator() );
        }
    }
}


