#ifndef _BUCKETFORSPHERES_H
#define _BUCKETFORSPHERES_H

#include "rg_Const.h"
#include "Sphere.h"
#include "rg_dList.h"
#include "AxisAlignedBox.h"
#include <math.h>

#include <set>
using namespace std;

template <class T>
class SphereElement
{
private:
    Sphere  m_sphere;
    T       m_property;
    rg_BOOL m_visited;

public:
    SphereElement();
    SphereElement(const Sphere& sphere, const T& sproperty);
    SphereElement(const SphereElement<T>& element);
    ~SphereElement();

    inline Sphere  getSphere() const { return m_sphere; }
    inline void    setSphere(const Sphere& sphere) { m_sphere = sphere; }

    inline T       getProperty() const { return m_property; }
    inline void    setProperty(const T& sproperty)  { m_property = sproperty; }

    inline rg_BOOL isVisited() const { return m_visited; }
    inline void    isVisited(const rg_BOOL& visited) { m_visited=visited; }

    SphereElement& operator =(const SphereElement<T>& element);
};


template <class T>
class BucketForSpheres
{
private:
    rg_dList< SphereElement<T> >     m_data;

    rg_dList< SphereElement<T>* >*** m_elements;
    rg_REAL                          m_min[3];
    rg_REAL                          m_max[3];
    rg_REAL                          m_unitSize;
    rg_INT                           m_size[3];

public:
    BucketForSpheres();
    BucketForSpheres(const BucketForSpheres& bucket);
    ~BucketForSpheres();

    rg_REAL          getUnitSize() const;
    AxisAlignedBox   getBoundingBox() const;

    BucketForSpheres& operator =(const BucketForSpheres& bucket);

    void    clean();
    void    construct( const rg_dList< pair<Sphere, T> >& values ,
                    const rg_REAL&                     unitSize );
    void    getBucketIndex(const rg_Point3D& pt, rg_INT& indexX, rg_INT& indexY, rg_INT& indexZ) const;


    void    pickOutSpheresByMask(   const rg_Point3D& pointMask, rg_dList<Sphere>& sphereList) const;
    void    pickOutPropertiesByMask(const rg_Point3D& pointMask, rg_dList<T>&      propertyList);

    void    pickOutSpheresByMask(   const Sphere& mask, rg_dList<Sphere>& sphereList) const;
    void    pickOutPropertiesByMask(const Sphere& mask, rg_dList<T>&      propertyList);

    void    pick_out_spheres_by_spherical_mask(const Sphere& mask, rg_dList< pair<Sphere, T> >& result) const;
    void    pick_out_spheres_to_contain_mask(const Sphere& mask, rg_dList< pair<Sphere, T> >& result, const rg_REAL& res = rg_SYSTEM_RES) const;


private:
    void    makeData( const rg_dList< pair<Sphere, T> >& values );
    void    makeBoundingBox( const rg_dList< pair<Sphere, T> >& values );
};



template<class T>
SphereElement<T>::SphereElement()
{
    m_visited = rg_FALSE;
}



template<class T>
SphereElement<T>::SphereElement(const Sphere& sphere, const T& sproperty)
{
    m_sphere   = sphere;
    m_property = sproperty;
    m_visited  = rg_FALSE;
}



template<class T>
SphereElement<T>::SphereElement(const SphereElement<T>& element)
{
    m_sphere   = element.m_sphere;
    m_property = element.m_property;
    m_visited  = element.m_visited;
}



template<class T>
SphereElement<T>::~SphereElement()
{
}



template<class T>
SphereElement<T>& SphereElement<T>::operator =(const SphereElement<T>& element)
{
    if ( this == &element ) {
        return *this;
    }

    m_sphere   = element.m_sphere;
    m_property = element.m_property;
    m_visited  = element.m_visited;

    return *this;
}



template<class T>
BucketForSpheres<T>::BucketForSpheres()
: m_elements(rg_NULL), m_unitSize(0.0)
{
    for ( rg_INT i=0; i<3; i++ ) {
        m_min[i]  = 0.0;
        m_max[i]  = 0.0;
        m_size[i] = 0;
    }
}



template<class T>
BucketForSpheres<T>::BucketForSpheres(const BucketForSpheres& bucket)
: m_unitSize(bucket.m_unitSize)
{
    rg_INT i = 0;
    rg_INT j = 0;
    rg_INT k = 0;
    for ( i=0; i<3; i++ ) {
        m_min[i]  = bucket.m_min[i];
        m_max[i]  = bucket.m_max[i];
        m_size[i] = bucket.m_size[i];
    }

    m_data.duplicateList( bucket.m_data );

    for ( i=0; i<m_size[0]; i++ ) {
        for ( j=0; j<m_size[1]; j++ ) {
            for ( k=0; k<m_size[2]; k++ ) {
                m_elements[i][j][k].duplicateList( bucket.m_elements[i][j][k] );
            }
        }
    }
}



template<class T>
BucketForSpheres<T>::~BucketForSpheres()
{
    clean();
}




template<class T>
rg_REAL BucketForSpheres<T>::getUnitSize() const
{
    return m_unitSize;
}



template<class T>
AxisAlignedBox BucketForSpheres<T>::getBoundingBox() const
{
    AxisAlignedBox boundingBox;
    boundingBox.setAxisAlignedBox( rg_Point3D(m_min[0], m_min[1], m_min[2]), rg_Point3D(m_max[0], m_max[1], m_max[2]) );

    return boundingBox;
}




template<class T>
BucketForSpheres<T>& BucketForSpheres<T>::operator =(const BucketForSpheres& bucket)
{
    if ( this == &bucket ) {
        return *this;
    }

    clean();

    m_unitSize   = bucket.m_unitSize;

    rg_INT i = 0;
    rg_INT j = 0;
    rg_INT k = 0;
    for ( i=0; i<3; i++ ) {
        m_min[i]  = bucket.m_min[i];
        m_max[i]  = bucket.m_max[i];
        m_size[i] = bucket.m_size[i];
    }

    m_data.duplicateList( bucket.m_data );

    for ( i=0; i<m_size[0]; i++ ) {
        for ( j=0; j<m_size[1]; j++ ) {
            for ( k=0; k<m_size[2]; k++ ) {
                m_elements[i][j][k].duplicateList( bucket.m_elements[i][j][k] );
            }
        }
    }
}




template<class T>
void BucketForSpheres<T>::clean()
{
    m_unitSize = 0.0;

    rg_INT i = 0;
    rg_INT j = 0;
    if ( m_elements != rg_NULL ) {
        for ( i=0; i<m_size[0]; i++ ) {
            for ( j=0; j<m_size[1]; j++ ) {
                delete [] m_elements[i][j];
            }
            delete [] m_elements[i];
        }
        delete [] m_elements;
    }

    for ( i=0; i<3; i++ ) {
        m_min[i]  = 0.0;
        m_max[i]  = 0.0;
        m_size[i] = 0;
    }

    m_data.removeAll();
}



template<class T>
void BucketForSpheres<T>::construct( const rg_dList< pair<Sphere, T> >& values,
                                     const rg_REAL&               unitSize )
{
    clean();


    makeData( values );
    makeBoundingBox( values );
    if (rg_NZERO(unitSize)) {
        m_unitSize = unitSize;
    }
    else {
        double totalVolume    = (m_max[0] - m_min[0])*(m_max[1] - m_min[1])*(m_max[2] - m_min[2]);
        double volumePerAData = totalVolume / m_data.getSize();
        m_unitSize = pow(volumePerAData, 1.0/3.0);
    }


	rg_INT i=0, j=0;

    for ( i=0; i<3; i++)  {
        rg_REAL size        = m_max[i] - m_min[i];
        rg_INT  numElements = (rg_INT) (size/m_unitSize);
        rg_REAL remnant     = size - (numElements*m_unitSize);
        rg_REAL increment   = m_unitSize - remnant;
        rg_REAL variation   = increment/2.;

        if ( remnant != 0.0 )  {
            m_max[i] += variation;
            m_min[i] -= variation;
        }

        rg_REAL numSize     = size/m_unitSize;
        rg_REAL ceilNumSize = ceil( numSize );

        if ( numSize < ceilNumSize && rg_ZERO(numSize - ceilNumSize) ) {
            m_size[i] = (rg_INT)( ceilNumSize );
        }
        else  {
            m_size[i] = (rg_INT)( numSize );
        }
    }

	m_elements = new rg_dList< SphereElement<T>* >**[m_size[0]];
	for( i=0; i<m_size[0]; i++)  {
		m_elements[i] = new rg_dList< SphereElement<T>* >* [m_size[1]];

        for( j=0; j<m_size[1]; j++)  {
			m_elements[i][j] = new rg_dList< SphereElement<T>* > [m_size[2]];
		}
	}

    //insert cell into bucket structure
    rg_Point3D center;
    rg_REAL    radius;

	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    m_data.reset4Loop();
    while ( m_data.setNext4Loop() ) {
        SphereElement<T>* ptrData = m_data.getpEntity();

        center = ptrData->getSphere().getCenter();
        radius = ptrData->getSphere().getRadius();

	    rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	    rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	    getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	    getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);

	    for(int i = minIndexX; i <=maxIndexX; i++) {
		    for(int j = minIndexY; j <= maxIndexY; j++) {
			    for(int k= minIndexZ; k <= maxIndexZ; k++) {
                    m_elements[i][j][k].add( ptrData );
			    }
		    }
	    }
	}
}



template<class T>
void BucketForSpheres<T>::makeData( const rg_dList< pair<Sphere, T> >& values )
{
    values.reset4Loop();
    while ( values.setNext4Loop() ) {
        pair<Sphere, T> data = values.getEntity();

        m_data.add( SphereElement<T>(data.first, data.second) );
    }
}



template<class T>
void BucketForSpheres<T>::makeBoundingBox( const rg_dList< pair<Sphere, T> >& values ) 
{
    for (rg_INT  i=0; i<3; i++ ) {
        m_min[i]  = DBL_MAX;
        m_max[i]  = -DBL_MAX;
    }

    rg_Point3D center;
    rg_REAL  x, y, z, radius;

    values.reset4Loop();
    while ( values.setNext4Loop() ) {
        pair<Sphere, T> data = values.getEntity();

        center = data.first.getCenter();
        x      = center.getX();
        y      = center.getY();
        z      = center.getZ();
        radius = data.first.getRadius();

        if ( (x + radius) > m_max[0] ) {
            m_max[0] = x + radius;
        }
        if ( (y + radius) > m_max[1] ) {
            m_max[1] = y + radius;
        }
        if ( (z + radius) > m_max[2] ) {
            m_max[2] = z + radius;
        }

        if ( (x - radius) < m_min[0] ) {
            m_min[0] = x - radius;
        }
        if ( (y - radius) < m_min[1] ) {
            m_min[1] = y - radius;
        }
        if ( (z - radius) < m_min[2] ) {
            m_min[2] = z - radius;
        }
    }
}



template<class T>
void    BucketForSpheres<T>::pickOutSpheresByMask(const rg_Point3D& pointMask, rg_dList<Sphere>& sphereList) const
{
    rg_INT indexX = -1;
    rg_INT indexY = -1;
    rg_INT indexZ = -1;
    getBucketIndex(pointMask, indexX, indexY, indexZ);


    m_elements[indexX][indexY][indexZ].reset4Loop();

    while (m_elements[indexX][indexY][indexZ].setNext4Loop()) {
        SphereElement<T>* element = m_elements[indexX][indexY][indexZ].getEntity();

        Sphere  sphere = element->getSphere();
        if ( sphere.doesContain(pointMask) ) {
            sphereList.add(sphere);
        }
    }
}



template<class T>
void    BucketForSpheres<T>::pickOutPropertiesByMask(const rg_Point3D& pointMask, rg_dList<T>&      propertyList)
{
    rg_INT indexX = -1;
    rg_INT indexY = -1;
    rg_INT indexZ = -1;
    getBucketIndex(pointMask, indexX, indexY, indexZ);


    m_elements[indexX][indexY][indexZ].reset4Loop();

    while (m_elements[indexX][indexY][indexZ].setNext4Loop()) {
        SphereElement<T>* element = m_elements[indexX][indexY][indexZ].getEntity();

        Sphere  sphere = element->getSphere();
        if (sphere.doesContain(pointMask)) {
            propertyList.add(element->getProperty());
        }
    }
}



template<class T>
void BucketForSpheres<T>::pickOutSpheresByMask(     const Sphere& mask, rg_dList<Sphere>& sphereList) const
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    rg_Point3D mask_center = mask.getCenter();
    rg_REAL    mask_radius = mask.getRadius();

	rg_Point3D sp = mask_center + rg_Point3D(-mask_radius, -mask_radius, -mask_radius);
	rg_Point3D ep = mask_center + rg_Point3D( mask_radius,  mask_radius,  mask_radius);

	getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);


    rg_dList< SphereElement<T>* > elementsPickedOut;

    rg_INT i, j, k;
	for( i = minIndexX; i <=maxIndexX; i++)  {
		for( j = minIndexY; j <= maxIndexY; j++)  {
			for( k= minIndexZ; k <= maxIndexZ; k++)  {

				m_elements[i][j][k].reset4Loop();
				
                while( m_elements[i][j][k].setNext4Loop() )  {
					SphereElement<T>* element = m_elements[i][j][k].getEntity();
                    
                    Sphere  sphere = element->getSphere();
                    rg_REAL distBtwCenters = mask_center.distance(sphere.getCenter());
                    rg_REAL sumOfTwoRadii  = mask_radius + sphere.getRadius();

					if( distBtwCenters < sumOfTwoRadii )  {
                        if ( !element->isVisited() ) {
                            element->isVisited( rg_TRUE );
                            elementsPickedOut.add(element);
                        }
                    }
				}
			}
		}
	}


    elementsPickedOut.reset4Loop();
    while ( elementsPickedOut.setNext4Loop() ) {
        SphereElement<T>* element = elementsPickedOut.getEntity();
        sphereList.add( element->getSphere() );
        element->isVisited(rg_FALSE);
    }
}



template<class T>
void BucketForSpheres<T>::pickOutPropertiesByMask( const Sphere& mask, rg_dList<T>& propertyList) 
{
	rg_INT maxIndexX = -1;
	rg_INT minIndexX = -1;
	rg_INT maxIndexY = -1;
	rg_INT minIndexY = -1;
	rg_INT maxIndexZ = -1;
	rg_INT minIndexZ = -1;

    rg_Point3D center = mask.getCenter();
    rg_REAL    radius = mask.getRadius();

	rg_Point3D sp = center + rg_Point3D(-radius, -radius, -radius);
	rg_Point3D ep = center + rg_Point3D( radius,  radius,  radius);

	getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
	getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);


    rg_dList< SphereElement<T>* > elementsPickedOut;
    rg_INT i, j, k;
	for( i = minIndexX; i <=maxIndexX; i++)  {
		for( j = minIndexY; j <= maxIndexY; j++)  {
			for( k= minIndexZ; k <= maxIndexZ; k++)  {

				m_elements[i][j][k].reset4Loop();
				
                while( m_elements[i][j][k].setNext4Loop() )  {
					SphereElement<T>* element = m_elements[i][j][k].getEntity();
                    
                    Sphere  sphere = element->getSphere();
                    rg_REAL distBtwCenters = center.distance(sphere.getCenter());
                    rg_REAL sumOfTwoRadii  = radius + sphere.getRadius();

					if( distBtwCenters < sumOfTwoRadii )  {
                        if ( !element->isVisited() ) {
                            element->isVisited( rg_TRUE );
                            elementsPickedOut.add(element);
                        }
                    }
				}
			}
		}
	}

    elementsPickedOut.reset4Loop();
    while ( elementsPickedOut.setNext4Loop() ) {
        SphereElement<T>* element = elementsPickedOut.getEntity();
        propertyList.add( element->getProperty() );
        element->isVisited(rg_FALSE);
    }
}



template<class T>
void BucketForSpheres<T>::pick_out_spheres_by_spherical_mask(const Sphere& mask, rg_dList< pair<Sphere, T> >& result) const
{
    rg_INT maxIndexX = -1;
    rg_INT minIndexX = -1;
    rg_INT maxIndexY = -1;
    rg_INT minIndexY = -1;
    rg_INT maxIndexZ = -1;
    rg_INT minIndexZ = -1;

    rg_Point3D mask_center = mask.getCenter();
    rg_REAL    mask_radius = mask.getRadius();

    rg_Point3D sp = mask_center + rg_Point3D(-mask_radius, -mask_radius, -mask_radius);
    rg_Point3D ep = mask_center + rg_Point3D(mask_radius, mask_radius, mask_radius);

    getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
    getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);


    rg_dList< SphereElement<T>* > elementsPickedOut;

    rg_INT i, j, k;
    for (i = minIndexX; i <= maxIndexX; i++) {
        for (j = minIndexY; j <= maxIndexY; j++) {
            for (k = minIndexZ; k <= maxIndexZ; k++) {

                m_elements[i][j][k].reset4Loop();

                while (m_elements[i][j][k].setNext4Loop()) {
                    SphereElement<T>* element = m_elements[i][j][k].getEntity();

                    Sphere  sphere = element->getSphere();
                    rg_REAL distBtwCenters = mask_center.distance(sphere.getCenter());
                    rg_REAL sumOfTwoRadii = mask_radius + sphere.getRadius();

                    if (distBtwCenters < sumOfTwoRadii) {
                        if (!element->isVisited()) {
                            element->isVisited(rg_TRUE);
                            elementsPickedOut.add(element);
                        }
                    }
                }
            }
        }
    }


    elementsPickedOut.reset4Loop();
    while (elementsPickedOut.setNext4Loop()) {
        SphereElement<T>* element = elementsPickedOut.getEntity();
        result.add( make_pair(element->getSphere(), element->getProperty()) );
        element->isVisited(rg_FALSE);
    }
}



template<class T>
void BucketForSpheres<T>::pick_out_spheres_to_contain_mask(const Sphere& mask, rg_dList< pair<Sphere, T> >& result, const rg_REAL& res) const
{
    rg_INT maxIndexX = -1;
    rg_INT minIndexX = -1;
    rg_INT maxIndexY = -1;
    rg_INT minIndexY = -1;
    rg_INT maxIndexZ = -1;
    rg_INT minIndexZ = -1;

    rg_Point3D mask_center = mask.getCenter();
    rg_REAL    mask_radius = mask.getRadius();

    rg_Point3D sp = mask_center + rg_Point3D(-mask_radius, -mask_radius, -mask_radius);
    rg_Point3D ep = mask_center + rg_Point3D(mask_radius, mask_radius, mask_radius);

    getBucketIndex(sp, minIndexX, minIndexY, minIndexZ);
    getBucketIndex(ep, maxIndexX, maxIndexY, maxIndexZ);


    rg_dList< SphereElement<T>* > elementsPickedOut;

    rg_INT i, j, k;
    for (i = minIndexX; i <= maxIndexX; i++) {
        for (j = minIndexY; j <= maxIndexY; j++) {
            for (k = minIndexZ; k <= maxIndexZ; k++) {

                m_elements[i][j][k].reset4Loop();

                while (m_elements[i][j][k].setNext4Loop()) {
                    SphereElement<T>* element = m_elements[i][j][k].getEntity();

                    Sphere  sphere = element->getSphere();
                    rg_REAL element_radius = sphere.getRadius();
                    rg_REAL distBtwCenters = mask_center.distance(sphere.getCenter());

                    if ( rg_LE( distBtwCenters + mask_radius, element_radius, res) ) {
                        if (!element->isVisited()) {
                            element->isVisited(rg_TRUE);
                            elementsPickedOut.add(element);
                        }
                    }
                }
            }
        }
    }


    elementsPickedOut.reset4Loop();
    while (elementsPickedOut.setNext4Loop()) {
        SphereElement<T>* element = elementsPickedOut.getEntity();
        result.add(make_pair(element->getSphere(), element->getProperty()));
        element->isVisited(rg_FALSE);
    }
}



template<class T>
void BucketForSpheres<T>::getBucketIndex(const rg_Point3D& pt, rg_INT& indexX, rg_INT& indexY, rg_INT& indexZ) const
{
	indexX = (rg_INT) ( ( pt.getX() - m_min[0] ) / ( m_max[0] - m_min[0] ) * m_size[0] );
	indexY = (rg_INT) ( ( pt.getY() - m_min[1] ) / ( m_max[1] - m_min[1] ) * m_size[1] );
	indexZ = (rg_INT) ( ( pt.getZ() - m_min[2] ) / ( m_max[2] - m_min[2] ) * m_size[2] );

	if( indexX > m_size[0]-1 ) indexX = m_size[0]-1;
	if( indexY > m_size[1]-1 ) indexY = m_size[1]-1;
	if( indexZ > m_size[2]-1 ) indexZ = m_size[2]-1;

	if( indexX < 0 ) indexX = 0;
	if( indexY < 0 ) indexY = 0;
	if( indexZ < 0 ) indexZ = 0;
}




#endif


