#ifndef _RG_NUBSPLINEINTERSECTOR_H
#define _RG_NUBSPLINEINTERSECTOR_H

#include "rg_Const.h"
#include "rg_dList.h"

#include "rg_NUBSplineCurve3D.h"

class rg_NUBSplineIntersector
{


public:
	rg_NUBSplineIntersector();
	~rg_NUBSplineIntersector();

	rg_dList<rg_Point2D> intersectBSplineCurveVsBSplineCurveUsingIntsecInterval(rg_NUBSplineCurve3D& curve_s, 
													                      rg_NUBSplineCurve3D& curve_t,
													                      rg_REAL & time);

	rg_dList<rg_Point2D> intersectBSplineCurveVsBSplineCurveUsingHardCodedImplicitization(rg_NUBSplineCurve3D& curve_s, 
													                                rg_NUBSplineCurve3D& curve_t,
													                                rg_REAL & time);

	rg_dList<rg_Point2D> intersectBSplineCurveVsBSplineCurveUsingImplicitization(rg_NUBSplineCurve3D& curve_s, 
													                       rg_NUBSplineCurve3D& curve_t,
													                       rg_REAL & time);
};

#endif



