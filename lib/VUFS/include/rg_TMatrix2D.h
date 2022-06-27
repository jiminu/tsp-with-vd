#ifndef _RG_TMATRIX2D_H
#define _RG_TMATRIX2D_H

#include "rg_Const.h"
#include "rg_Matrix.h"
#include "rg_Point2D.h"
#include "rg_Line2D.h"

class rg_TMatrix2D :public rg_Matrix
{
public:
	rg_TMatrix2D();
	rg_TMatrix2D(const rg_INT &dg);
	rg_TMatrix2D(const rg_TMatrix2D &matrix);
	~rg_TMatrix2D();

	void scale(const rg_REAL &a, const rg_REAL &b);
	void translate(const rg_REAL &m, const rg_REAL &n);
    void translate(const rg_Point2D& pt);
//	void rotate(const rg_REAL &m, const rg_REAL &n, const rg_REAL &theta);
    void rotate(const rg_REAL& angle);
    void rotate(const rg_REAL& angle, const rg_Point2D& center);
    void rotate(const rg_Point2D& fromVector, const rg_Point2D& endVector);
    void rotate(const rg_Point2D& toVector);
//	void origin_rotate(const rg_REAL &theta);

    void reflectToXAxis();
    void reflectToYAxis();
    void reflect(const rg_Line2D& line);

	rg_TMatrix2D operator*(const rg_TMatrix2D &matrix) const;
	rg_Point2D operator*(const rg_Point2D &point) const;
	rg_TMatrix2D& operator=(const rg_TMatrix2D &matrix);
};

///////////////////////////
#endif


