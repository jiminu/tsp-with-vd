/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_Matrix.h
//	  
//    DESCRIPTION : 
//           This is the interface of the class rg_Matrix 
//           which define rg_Matrix and its property. 
//                          
//	  CLASS NAME  : rg_Matrix
//
//    BASE CLASS  : None
//      
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 1 Jul. 1997    
//    
//    HISTORY :
//
//           1. Lee, Soon-Woong inserted the following function at 19 Jul. 1997.
//	             void     gaussianElimination()
//	             void     fwdElimination()
//	             void     bkSubstitution()
//	             rg_REAL   determinant()
//	             void     LUdecomposition()
//
//           2. Jang, Tae-Bum inserted the following function at 1997.10.3 due to OPEN_GL
//               rg_REAL*  getColMajorOrderedArray() const
//               rg_REAL*  getRowMajorOrderedArray() const
//
//           3. Ryu, Jung-Hyun 18. Jul. 1998.
//               rg_Matrix getOneRowDeletedSubmatrix(const rg_INDEX& theCorrespondingRow) const;
//               rg_Matrix getOneColDeletedSubmatrix(const rg_INDEX& theCorrespondingCol) const;
//
//           4. Cho, Youngsong 23. Apr. 1999.
//			     rg_Matrix inverse()  const;   //  Gauss-Jordan elimination
//               rg_FLAG   isIdentity() const;
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _MATRIX_H
#define _MATRIX_H

#include <math.h>
#include <string.h>
#include "rg_Const.h"
//#include <fstream>


class rg_Matrix
{
protected:
	rg_REAL **element;
	rg_INT row;
	rg_INT col;

public:
	rg_Matrix();
	rg_Matrix(const rg_INT &r, const rg_INT &c);
	rg_Matrix(const rg_Matrix& matrix);
	~rg_Matrix();

	void removeAll();
	void setSize(const rg_INT &r,const rg_INT &c);
    void setElement(const rg_INT &r, const  rg_INT &c,const rg_REAL &value);
	void setElements(rg_REAL** tElement, const rg_INT& rowSize, const rg_INT& colSize);
	void makeIdentity();
	rg_FLAG isIdentity() const;

	rg_REAL** getElement() const;
	rg_REAL getElement(const rg_INT &r,const rg_INT &c) const;
	rg_INT getColSize()const ;
	rg_INT getRowSize()const ;

	rg_REAL trace() const;
	rg_Matrix transpose() const;
	//  Gauss-Jordan elimination
	rg_Matrix inverse() const;


	rg_REAL  getNorm(const rg_INT &k) const;
	rg_REAL  getNorm(const char* ch) const;
	rg_Matrix getOneRowDeletedSubmatrix(const rg_INDEX & theCorrespondingRow) const;
	rg_Matrix getOneColDeletedSubmatrix(const rg_INDEX & theCorrespondingCol) const;

	rg_REAL *operator[](const rg_INT &r);
	rg_Matrix operator*(const rg_Matrix &matrix)const;
	rg_Matrix operator*(const rg_REAL &n) const;
	rg_Matrix operator/(const rg_REAL &n) const;
	rg_Matrix operator+(const rg_Matrix &matrix) const;
	rg_Matrix operator-(const rg_Matrix &matrix) const;
	rg_Matrix& operator=(const rg_Matrix &matrix) ;

//	void show();

	friend rg_Matrix operator*(const rg_REAL &scalar,const rg_Matrix& matrix);

// The following was inserted by Lee, Soon-Woong at 19 Jul. 1997.
	void gaussianElimination(rg_Matrix &x, rg_Matrix &b);
	void fwdElimination(rg_Matrix &m, rg_Matrix &b) const;
	void bkSubstitution(rg_Matrix &m, rg_Matrix &x, rg_Matrix &b) const;
	rg_REAL determinant() const;
	void LUdecomposition(rg_Matrix &L, rg_Matrix &U);
// End of Insertion       

// The following was inserted by Jang, Tae-Bum at 1997.10.30 for OPEN_GL
//    ColMajorOrdered is abbrieviation of "Column-major ordered".
//    This funtion retun following array
//    | a00 a01 a02 a03 |
//    | a10 a11 a12 a13 |
//    | a20 a21 a22 a23 |   -> [ a00, a10, a20, a30, a10, a11,......,a33]
//    | a30 a31 a32 a33 |
//
//    In OPEN_GL matrix is stored above.
//     * Attention ; must deleted return pointer, which is allocated by new
    rg_REAL*  getColMajorOrderedArray() const;
//
//    In contrast, 
//    "Row major orderd" storage is to store the elements of matrix as follows
//    ->[ a00, a01, a02, a03,..........., a33]
    rg_REAL*  getRowMajorOrderedArray() const;
// End of Insertion
	void display();
/*
	friend ofstream& operator<<(ofstream& fout, const rg_Matrix& mat);
    friend std::ostream& operator<<(std::ostream& out, const rg_Matrix& mat);
*/
};

#endif
