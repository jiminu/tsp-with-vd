/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_EllipticArc3D.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_EllipticArc3D 
//           which define ellipse and its property. 
//
//    BASE CLASS  : None  
//
//    AUTHOR      : Taeboom Jang
//    START DATE  : 2001.  1.  29. 
//
//    HISTORY     :
//
//    
//            Copyright (c) CAD/CAM Lab.    	  	
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_ELLIIPTICARC3D_H_
#define _RG_ELLIIPTICARC3D_H_

#include "rg_Const.h"
#include "rg_Point3D.h"

class rg_EllipticArc3D
{
private:
	rg_Point3D centerPt;
    rg_Point3D localXAxis;
    rg_Point3D localYAxis;
    rg_REAL    startAngle;
    rg_REAL    endAngle;
	
public:
	rg_EllipticArc3D();
	rg_EllipticArc3D(const rg_EllipticArc3D &temp);
	~rg_EllipticArc3D();
	
//get functions
	rg_Point3D  getCenterPt()   const;
	rg_Point3D  getLocalXAxis() const;
	rg_Point3D  getLocalYAxis() const;
	rg_REAL     getStartAngle() const;
	rg_REAL     getEndAngle()   const;
	
//set functions
	void setCenterPt(const rg_Point3D &center);
	void setLocalXAxis(const rg_Point3D &tLocalXAxis);
	void setLocalYAxis(const rg_Point3D &tLocalYAxis);
	void setStartAngle(const rg_REAL &tStartAngle);
	void setEndAngle(const rg_REAL &tEndAngle);
    void setEllipticArc3D(const rg_EllipticArc3D &temp);

//evaluate
    rg_Point3D evaluatePt(const rg_REAL& angle) const;
    rg_Point3D evaluateDerivative(const rg_REAL& angle) const;

// ETC
    rg_FLAG     isOn(const rg_Point3D& pt) const;
    rg_Point3D* getFoci() const;
    rg_Point3D  getMajorAxis() const;
    rg_Point3D  getMinorAxis() const;
    rg_REAL     getAngle(const rg_Point3D& pt) const;
//operators	
	rg_EllipticArc3D& operator=(const rg_EllipticArc3D& temp);

};

#endif


