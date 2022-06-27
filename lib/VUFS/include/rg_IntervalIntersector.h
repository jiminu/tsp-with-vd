
#ifndef rg_INTERVALINTERSECTOR_H
#define rg_INTERVALINTERSECTOR_H

#include "rg_Const.h"
#include "rg_BzCurve2D.h"
#include "rg_Box2DWithParameter.h"
#include "rg_dList.h"
#include "rg_Polynomial.h"
class rg_Point2D;

class rg_IntervalIntersector
{
private:
	rg_BzCurve2D            curve_s;
	rg_BzCurve2D            curve_t;
	rg_dList<rg_Point2D> intersectionPointList;

public:
	rg_IntervalIntersector();
	rg_IntervalIntersector(const rg_BzCurve2D &curve1, const rg_BzCurve2D &curve2, const rg_dList<rg_Point2D> &ptlist);
	~rg_IntervalIntersector();

	void		   intersectBzCurveVsBzCurve(const rg_BzCurve2D& curve1, const rg_BzCurve2D& curve2, rg_REAL &time);
	rg_dList<rg_Box2DWithParameter>     predetermineBox(const rg_BzCurve2D& curve);
	void           boxingLoop(rg_Box2DWithParameter box_s, rg_Box2DWithParameter box_t);
	rg_INT           isOverlapping(rg_Box2DWithParameter box_s, rg_Box2DWithParameter box_t);
	void           changingBox(rg_Box2DWithParameter &box);
	rg_INT           isTerminalCondition(const rg_Box2DWithParameter &box_s, 
								    const rg_Box2DWithParameter &box_t, 
									rg_Point2D &point);
	rg_Box2DWithParameter*           splitBox(const rg_Box2DWithParameter& box, const rg_BzCurve2D& curve);
	rg_dList<rg_Point2D> getIntersectionPointList();
	rg_dList<rg_REAL> getLocalMaxMinPointsOntheXY(const rg_BzCurve2D& curve);
//	rg_Polynomial convertBzCurve2Polynomial(const rg_BzCurve2D& curve, const rg_REAL* b);
};

#endif // rg_IntervalIntersector class


