#include "Generator2D.h"
#include "VFace2D.h"
#include "VoronoiDiagramCIC.h"
using namespace V::GeometryTier;


Generator2D::Generator2D()
    : m_type(Generator_Type::DISK_G)
    , m_ID(0)
    , m_outerFace(NULL)
    , m_innerFace(NULL)
    , m_innerVD(NULL)
    , m_userData(NULL)
{
}



Generator2D::Generator2D( const int& ID )
    : m_type(Generator_Type::DISK_G)
    , m_ID(ID)
    , m_outerFace(NULL)
    , m_innerFace(NULL)
    , m_innerVD(NULL)
    , m_userData(NULL)
{
}



Generator2D::Generator2D( const int& ID, const rg_Circle2D& disk )
    : m_type(Generator_Type::DISK_G)
    , m_ID(ID)
    , m_disk(disk)
    , m_outerFace(NULL)
    , m_innerFace(NULL)
    , m_innerVD(NULL)
    , m_userData(NULL)
{
}



Generator2D::Generator2D( const int& ID, const rg_Circle2D& disk, VFace2D* const face )
    : m_type(Generator_Type::DISK_G)
    , m_ID(ID)
    , m_disk(disk)
    , m_outerFace(face)
    , m_innerFace(NULL)
    , m_innerVD(NULL)
    , m_userData(NULL)
{
}



Generator2D::Generator2D( const int& ID, const rg_Circle2D& disk, void* const userData )
    : m_type(Generator_Type::DISK_G)
    , m_ID(ID)
    , m_disk(disk)
    , m_outerFace(NULL)
    , m_innerFace(NULL)
    , m_innerVD(NULL)
    , m_userData(userData)
{
}



Generator2D::Generator2D( const Generator2D& generator )
    : m_type(       generator.m_type)
    , m_ID(         generator.m_ID)
    , m_disk(       generator.m_disk)
    , m_outerFace(  generator.m_outerFace)
    , m_innerFace(  generator.m_innerFace)
    , m_innerVD(    generator.m_innerVD)
    , m_userData(   generator.m_userData)
{

    // CYSONG ADDED THE FOLLOWINGS. (SEP 16, 20)
    m_innerGens = generator.m_innerGens;
}



Generator2D::Generator2D( const Generator_Type& type )
    : m_type(type)
    , m_ID(0)
    , m_outerFace(NULL)
    , m_innerFace(NULL)
    , m_innerVD(NULL)
    , m_userData(NULL)
{
}



Generator2D::Generator2D( const Generator_Type& type, const int& ID )
    : m_type(type)
    , m_ID(ID)
    , m_outerFace(NULL)
    , m_innerFace(NULL)
    , m_innerVD(NULL)
    , m_userData(NULL)
{
}



Generator2D::Generator2D( const Generator_Type & type, const VFace2D * const Vface, void * userData, const int& ID )
    : m_type(type)
    , m_ID(ID)
    , m_outerFace(const_cast<VFace2D*>(Vface))
    , m_innerFace(NULL)
    , m_innerVD(NULL)
    , m_userData(NULL)
{
}



Generator2D::~Generator2D()
{
//     if( m_innerVD != NULL )
//         delete m_innerVD;

    if (m_innerVD != NULL)
    {
        m_innerVD->setContainerGenerator(NULL);
        m_innerVD->clear_except_current_container(); 

        if (m_outerFace != NULL) 
        {
            delete m_innerVD;   
        }
    }
}



bool Generator2D::constructVD()
{
    if ( m_innerGens.size() < 2 ) {
        return false;
    }

    if ( m_innerVD != NULL ) 
        delete m_innerVD;

    m_innerVD = new VoronoiDiagramCIC();
    m_innerVD->constructVoronoiDiagramCIC(this);

    return true;
}



Generator2D& Generator2D::operator=( const Generator2D& generator )
{
    if( this == &generator ) {
        return *this;
    }

    m_type          = generator.m_type;
    m_ID            = generator.m_ID;
    m_disk          = generator.m_disk;
    m_outerFace     = generator.m_outerFace;
    m_innerFace     = generator.m_innerFace;
    m_innerVD       = generator.m_innerVD;
    m_userData      = generator.m_userData;

    return *this;
}




bool Generator2D::operator==( const Generator2D& generator ) const
{
    if( this == &generator ) {
        return true;
    }

    if(     m_type          == generator.m_type
        &&  m_ID            == generator.m_ID 
        &&  m_disk          == generator.m_disk 
        &&  m_outerFace     == generator.m_outerFace
        &&  m_innerFace     == generator.m_innerFace
        &&  m_innerVD       == generator.m_innerVD
        &&  m_userData      == generator.m_userData) 
    {
        return true;
    }
    else {
        return false;
    }
}



void Generator2D::getNeighborGeneratorsInThisVoronoiDiagram( list<Generator2D*>& neighborGeneratorsList, VoronoiDiagram2DC* const VD )
{
    list<VFace2D*> adjacentVFaces;
	if ( m_innerVD == VD )
        m_innerFace->getAdjacentVFaces( adjacentVFaces );
    else
        m_outerFace->getAdjacentVFaces( adjacentVFaces );

    for( list<VFace2D*>::iterator i_face = adjacentVFaces.begin() ; i_face != adjacentVFaces.end() ; ++i_face ) {
        VFace2D* currVFace = *i_face;

        if( currVFace->getGenerator() != NULL )
            neighborGeneratorsList.push_back( currVFace->getGenerator() );
    }
}



void Generator2D::getNeighborGenerators(list<Generator2D*>& neighborGeneratorsList, const bool & isThisContainer) const
{
    list<VFace2D*> adjacentVFaces;
    if (isThisContainer) {
        m_innerFace->getAdjacentVFaces(adjacentVFaces);
    }
    else {
        m_outerFace->getAdjacentVFaces(adjacentVFaces);
    }

    for (list<VFace2D*>::iterator i_face = adjacentVFaces.begin(); i_face != adjacentVFaces.end(); ++i_face) {
        VFace2D* currVFace = *i_face;

        if (currVFace->getGenerator() != NULL) {
            neighborGeneratorsList.push_back(currVFace->getGenerator());
        }
    }
}


