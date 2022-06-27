#ifndef _T3DVERTEX_H
#define _T3DVERTEX_H

#include "TopologicalEntity.h"
#include "rg_Point3D.h"

namespace V {

namespace GeometryTier {


const rg_INT NOT_VISITED = 0;
const rg_INT ON_PROCESS  = 2;
const rg_INT COMPLETED   = 1;

class T3DTetrahedron;

class T3DVertex : public TopologicalEntity
{
private:
	rg_Point3D      m_point;
	T3DTetrahedron* m_firstTetrahedron;
    rg_FLAG         m_check;
    void*           m_property;

public:
	T3DVertex();
	T3DVertex(const rg_Point3D& point);
	T3DVertex(const rg_INT& ID, const rg_Point3D& point);
	T3DVertex(const T3DVertex& vertex);
	~T3DVertex();

	rg_Point3D      getPoint() const;
	T3DTetrahedron* getFirstTetrahedron() const;    
    void*           getProperty() const;

    rg_FLAG isChecked() const;

	void setPoint(const rg_Point3D& point);
	void setFirstTetrahedron(T3DTetrahedron* firstTetrahedron);
    void setCheck(const rg_FLAG& check);
	void setProperty(void* property);
	void setVertex(const rg_Point3D& point, T3DTetrahedron* firstTetrahedron, void* property);

	T3DVertex& operator =(const T3DVertex& vertex);

};

} // namespace GeometryTier

} // namespace V

#endif 
