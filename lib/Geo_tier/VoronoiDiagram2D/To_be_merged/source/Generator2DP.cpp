#include "Generator2DP.h"
#include "VFace2DP.h"
using namespace BULL2D::GeometryTier;


Generator2DP::Generator2DP()
{
    m_ID            = -1;
    m_location      = NULL;
    m_assignedVFace = NULL;
}



Generator2DP::Generator2DP( const int& ID )
{
    m_ID            = ID;
    m_location      = NULL;
    m_assignedVFace = NULL;
}



Generator2DP::Generator2DP( const int& ID, rg_Point2D* const location )
{
    m_ID            = ID;
    m_location      = location;
    m_assignedVFace = NULL;
}



Generator2DP::Generator2DP( const int& ID, rg_Point2D* const location, VFace2DP* const face )
{
    m_ID            = ID;
    m_location      = location;
    m_assignedVFace = face;
}



Generator2DP::Generator2DP( const Generator2DP& generator )
{
    m_ID            = generator.m_ID;
    m_location      = generator.m_location;
    m_assignedVFace = generator.m_assignedVFace;
}



Generator2DP::~Generator2DP()
{
}



Generator2DP& Generator2DP::operator=( const Generator2DP& generator )
{
    if( this == &generator ) {
        return *this;
    }

    m_ID            = generator.m_ID;
    m_location      = generator.m_location;
    m_assignedVFace = generator.m_assignedVFace;

    return *this;
}



void Generator2DP::getNeighborGenerators( list<Generator2DP*>& neighborGeneratorsList )
{
    list<VFace2DP*> adjacentVFaces;
    m_assignedVFace->getAdjacentVFaces( adjacentVFaces );

    for( list<VFace2DP*>::iterator it_face = adjacentVFaces.begin() ; it_face != adjacentVFaces.end() ; it_face++ ) {
        VFace2DP* currVFace = *it_face;

        if( currVFace->getGenerator() != NULL ) {
            neighborGeneratorsList.push_back( currVFace->getGenerator() );
        }
    }
}


