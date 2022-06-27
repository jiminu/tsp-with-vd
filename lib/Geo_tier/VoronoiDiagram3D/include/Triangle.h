#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include "rg_Const.h"
#include "rg_Point3D.h"

namespace V {

namespace GeometryTier {


class Triangle
{
protected:
	rg_Point3D m_vertex[3];
	rg_Point3D m_normalOnVertex[3];

public:
	Triangle();
	Triangle(const Triangle& aTriangle);
	~Triangle();

	rg_Point3D getVertex(const rg_INT& i) const;
	rg_Point3D getNormalOnVertex(const rg_INT& i) const;

	void setVertex(const rg_INT& i, const rg_Point3D& aVertex);
	void setNormalOnVertex(const rg_INT& i, const rg_Point3D& aNormal);

	Triangle& operator =(const Triangle& aTriangle);
};


} // namespace GeometryTier

} // namespace V


#endif
