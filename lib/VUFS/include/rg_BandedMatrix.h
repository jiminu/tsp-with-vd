//	This is a part of the HANYANG CAD/CAM Lab. Library
//
//	Hanyang University 		
//	Department of Industrial Engineering		
//	Copyright ¨Ï 1996 by HANYANG CAD/CAM Lab.
//						
//	Written in Korea				
//								
//	All rights reserved. No part of this program may be reproduced 
//	or transmitted in any form or by any means without permission 
//	in programming from CAD/CAM Laboratory.		
//							
//	This source code was written by Lee SoonWoong in Aug. 17. 1996

#ifndef _BANDEDMATRIX_H_TEMPLATE
#define _BANDEDMATRIX_H_TEMPLATE	
#include "rg_Const.h"

template <class T>
class rg_BandedMatrix
{
private:
	rg_INT    matrixSize;
	rg_INT    bandSize;
	short  invertFlag;

	void   fwdEliminate(rg_REAL **m, T *b);
	void   bkSubstitute(rg_REAL **m, T *x, T *b);

public:
	rg_BandedMatrix();
	rg_BandedMatrix(const rg_INT &mSize, const rg_INT &bSize);
	~rg_BandedMatrix();

	void findVariable(rg_REAL **m, T *x, T *b);
	rg_REAL **convertMatrix(rg_REAL **tempm);
};

#endif

///////////////////////////////////////////////////////////
//  Implementation of template class

template <class T>
rg_BandedMatrix<T>::rg_BandedMatrix()
{
	matrixSize = 0;
	bandSize = 0;
	invertFlag = 0;
}

template <class T>
rg_BandedMatrix<T>::rg_BandedMatrix(const rg_INT &mSize, const rg_INT &bSize)
{
	matrixSize = mSize;
	bandSize = bSize;
	invertFlag = 0;
}

template <class T>
rg_BandedMatrix<T>::~rg_BandedMatrix()
{

}

template <class T>
void rg_BandedMatrix<T>::findVariable(rg_REAL **m, T *x, T *b)
{
	/*
	rg_REAL **m;
	if( sizeof(tempm[0]) == sizeof(tempm[1])) m=convertMatrix(tempm);
	else {
			m = new rg_REAL*[bandSize*2+1];
			for(rg_INT i=0; i<=bandSize; i++)
				m[i] = new rg_REAL[matrixSize-bandSize+i];
			for(i=bandSize+1; i<=2*bandSize; i++)
				m[i] = new rg_REAL[matrixSize+bandSize-i];

			for(i=0; i<=bandSize; i++)
				for(rg_INT j=0; j<matrixSize-bandSize+i; j++)
					m[i][j] = tempm[i][j];

			for(i=bandSize+1; i<=2*bandSize; i++)
				for(rg_INT j=0; j<matrixSize+bandSize-i; j++)
					m[i][j] = tempm[i][j];
	}
	*/

	fwdEliminate(m, b);
	bkSubstitute(m, x, b);

	if(invertFlag) 
	{
		rg_INT i = 0;
		for(i=0; i<bandSize*2+1; i++) delete [] m[i];
		delete [] m;
	}
}

template <class T>
rg_REAL** rg_BandedMatrix<T>::convertMatrix(rg_REAL **tempm)
{
	rg_REAL **m = new rg_REAL *[bandSize*2+1];
	invertFlag = 1;

	rg_INT i = 0;
	for(i=0; i<=bandSize; i++)
		m[i] = new rg_REAL[matrixSize-bandSize+i];
	for(i=bandSize+1; i<=2*bandSize; i++)
		m[i] = new rg_REAL[matrixSize+bandSize-i];

	for(i=0; i<=bandSize; i++)
		for(rg_INT j=0; j<matrixSize-bandSize+i; j++)
			m[i][j] = tempm[bandSize-i+j][j];

	for(i=bandSize+1; i<=2*bandSize; i++)
		for(rg_INT j=0; j<matrixSize+bandSize-i; j++)
			m[i][j] = tempm[j][i+j-bandSize];

	return m;
}

template <class T>
void rg_BandedMatrix<T>::fwdEliminate(rg_REAL **m, T *b)
{
	rg_REAL temp;
	rg_INT    tempi;

	for(rg_INT i=1, s=0; i<=matrixSize-1; i++) 
	{
		if(i>matrixSize-bandSize) s++;
		for(rg_INT j=0; j<=bandSize-s-1; j++) 
		{
			temp = m[bandSize-j-1][i-1]/m[bandSize][i-1];
			b[i+j] -= temp * b[i-1];
			for(rg_INT k=0; k<=bandSize-s-1; k++) 
			{
				if(j>=k) tempi = i+k;
				m[bandSize-j+k][tempi] -= temp * m[bandSize+k+1][i-1];
			}
		}
	}
}

template <class T>
void rg_BandedMatrix<T>::bkSubstitute(rg_REAL **m, T *x, T *b)
{
	T tempsum;
	rg_INT    s;

	for(rg_INT i=matrixSize-1; i>=0; i--) 
	{
		tempsum = 0;

		if(i>=matrixSize-bandSize) s=matrixSize-i;
		else s=bandSize+1;

		for(rg_INT j=1; j<s; j++)
			tempsum += m[bandSize+j][i]*x[i+j];
		x[i] = (b[i] - tempsum)/m[bandSize][i];
	}
}


