#include "MeshForSphere.h"
using namespace V::GeometryTier;

MeshForSphere::MeshForSphere()
: m_resolution(0)
{
}



MeshForSphere::MeshForSphere(const rg_INT& resolution, const Sphere& sphere)
: m_resolution(resolution), m_sphere(sphere)
{
    makeSphericalTriangularMesh(m_resolution, m_sphere.getCenter(), m_sphere.getRadius());
}



//MeshForSphere::MeshForSphere(const MeshForSphere& sphereMesh)
//{
//}



MeshForSphere::~MeshForSphere()
{

}




Sphere MeshForSphere::getSphere() const
{
    return m_sphere;
}



void   MeshForSphere::setSphere(const Sphere& sphere)
{
    translateSphericalTriangularMesh( -m_sphere.getCenter() );

    scaleSphericalTriangularMesh( sphere.getRadius() );
    translateSphericalTriangularMesh( sphere.getCenter() );

    m_sphere = sphere;
}



void   MeshForSphere::setResolution(const rg_INT& resolution)
{
    m_resolution = resolution;

    makeSphericalTriangularMesh(m_resolution, m_sphere.getCenter(), m_sphere.getRadius());
}




//MeshForSphere& MeshForSphere::operator =(const MeshForSphere& sphereMesh)
//{
//}




