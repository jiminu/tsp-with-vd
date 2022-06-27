/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_ImplicitEqMatrix.h
//	  
//    DESCRIPTION : 
//           This is the implementation of the class ImplicitMatrix
//           and is used for Implicitization method(Generating Implicit Eq.) alone!!
//
//    AUTHOR      : Ryu, JungHyun
//    START DATE  : 17 Feb. 1998    
//
//    HISTORY     :
//                 1.	 98. 3. 3 rg_ImplicitEqMatrix deleteOneElement(const rg_INDEX& theCorrespondingRow
//	                                                                const rg_INDEX& theCorrespondingCol) const;
//                                -> this member function replaces functions "deleteOneRow & deleteOneCol" (for the sake of efficiency)
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#ifndef _RG_EMPLICITEQMATRIX_IS_H
#define _RG_EMPLICITEQMATRIX_IS_H

#include "rg_Const.h"
#include "rg_ImplicitEquation.h"

class rg_ImplicitEqMatrix
{
public:
	rg_ImplicitEquation** element;
	rg_INT row;
	rg_INT col;

	// constructors and destructor
public:	
	rg_ImplicitEqMatrix();
	//rg_ImplicitEqMatrix(const rg_DEGREE tdegree);
	rg_ImplicitEqMatrix(const rg_INT &trow, const rg_INT &tcol);
	~rg_ImplicitEqMatrix();
	rg_ImplicitEqMatrix(const rg_ImplicitEqMatrix& SourceObj);

	// member functions and overloaded operators
public:
	rg_INT getRowSize() const;
	rg_INT getColSize() const;
	rg_ImplicitEquation getCofactor(const rg_INDEX& i, const rg_INDEX& j) const;
	rg_ImplicitEquation getDeterminant() const;
	rg_ImplicitEqMatrix deleteOneRow(const rg_INDEX& theCorrespondingRow) const;//antique
	rg_ImplicitEqMatrix deleteOneCol(const rg_INDEX& theCorrespondingCol) const;//antique

	rg_ImplicitEqMatrix deleteOneElement(const rg_INDEX& theCorrespondingRow,
	                                  const rg_INDEX& theCorrespondingCol) const;
	void             removeAll();
	rg_ImplicitEquation getElement(const rg_INDEX& i, const rg_INDEX& j) const;

	rg_ImplicitEquation* operator[](const rg_INDEX& index);
	rg_ImplicitEqMatrix& operator=(const rg_ImplicitEqMatrix& SourceObj);

	void show() const;
};

#endif


