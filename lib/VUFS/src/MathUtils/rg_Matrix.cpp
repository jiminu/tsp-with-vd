/////////////////////////////////////////////////////////////////////
//
//    FILENAME    : rg_Matrix.cpp
//    
//    DESCRIPTION : 
//           This is the implementation of the class rg_Matrix 
//
//    AUTHOR      : Lee, Soon-Woong
//    START DATE  : 8 Jun. 1997    
//
//    HISTORY     :
//
//           1. Lee, Soon-Woong inserted the following function at 19 Jul. 1997.
//                      rg_REAL determinant()
//               void     gaussianElimination()
//               void     fwdElimination()
//               void     bkSubstitution()
//               rg_REAL determinant()
//               void     LUdecomposition()
//			 2. Joonghyun Ryu 20 May. 2005
//				 rg_REAL trace()
//
//           Copyright ¨Ï 1996 by CAD/CAM Lab. in Hanyang University
//
/////////////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "rg_RelativeOp.h"

#include "rg_Matrix.h"

rg_Matrix::rg_Matrix()
{
    col     = 0;
    row     = 0;
    element = rg_NULL;
}

rg_Matrix::rg_Matrix(const rg_INT &r, const rg_INT &c)
{
    row = r;
    col = c;

    element = new rg_REAL*[row];
    rg_INT i = 0;
	for(i = 0; i < row; i++)
        element[i] = new rg_REAL[col];
    
    for(i = 0; i < row; i++)
        for(rg_INT j = 0; j < col; j++)
            element[i][j] = 0.0;
}

rg_Matrix::rg_Matrix(const rg_Matrix& matrix)
{
    row = matrix.row;
    col = matrix.col;

    element = new rg_REAL*[row];
    rg_INT i = 0;
	for(i = 0; i < row; i++)
        element[i] = new rg_REAL[col];
    
    for(i = 0; i < row; i++)
        for(rg_INT j = 0; j < col; j++)
            element[i][j] = matrix.element[i][j];
}

rg_Matrix::~rg_Matrix()
{
    removeAll();
}


void rg_Matrix::removeAll()
{
    if( element != rg_NULL ) 
    {
        for(rg_INT i=0; i<row; i++)
            delete element[i];

        delete []element;
    }

    element = rg_NULL;
    row     = 0;
    col     = 0;
}

void rg_Matrix::setSize(const rg_INT &r,const rg_INT &c)
{
    removeAll();

    row = r;
    col = c;

    element = new rg_REAL*[row];
    rg_INT i = 0;
	for(i = 0; i < row; i++)
        element[i] = new rg_REAL[col];
    
    for(i = 0; i < row; i++)
        for(rg_INT j = 0; j < col; j++)
            element[i][j] =0;
}

void rg_Matrix::setElement(const rg_INT &r, const  rg_INT &c,const rg_REAL &value)
{
   element[r][c] = value;
}

void rg_Matrix::setElements(rg_REAL** tElement, const rg_INT& rowSize, const rg_INT& colSize)
{
    setSize(rowSize, colSize);

    for (int i = 0; i < row; i++)
        for (rg_INT j = 0; j < col; j++)
            element[i][j] = tElement[i][j];
}

void rg_Matrix::makeIdentity()
{
    if(row != col) 
    {
        return;
    }

    for(rg_INT i = 0; i < row; i++)
    {
        for(rg_INT j = 0; j < col; j++)
        {
            if(i==j) element[i][j] = 1.0;
            else element[i][j] = 0.0;
        }
    }
}

rg_FLAG rg_Matrix::isIdentity() const
{
	if ( row != col )  {
		return rg_FALSE;
	}

	for (rg_INT i=0; i<row; i++)  {
		for (rg_INT j=0; j<col; j++)  {
			if ( i == j )  {
				if ( rg_NE(element[i][j], 1.0) )  {
					return rg_FALSE;
				}
			}
			else  {
				if ( rg_NZERO ( element[i][j] ) )  {
					return rg_FALSE;
				}
			}
		}
	}

	return rg_TRUE;
}

rg_REAL** rg_Matrix::getElement() const
{
	return element;
}

rg_REAL rg_Matrix::getElement(const rg_INT &r,const rg_INT &c) const
{
    return element[r][c];
}

rg_INT rg_Matrix::getColSize() const
{
    return col;
}

rg_INT rg_Matrix::getRowSize() const
{
    return row;
}

rg_REAL rg_Matrix::trace() const
{
	rg_REAL trace = 0.;

	for(rg_INT i = 0;i < row;i++)
		trace += element[ i ][ i ];

	return trace;
}

rg_Matrix rg_Matrix::transpose() const
{
    rg_Matrix transpz(col, row);
    
    for(rg_INT i=0; i<col; i++)
        for(rg_INT j=0; j<row; j++)
            transpz.element[i][j] = element[j][i];

    return transpz;
}

//  Gauss-Jordan elimination
rg_Matrix rg_Matrix::inverse()  const
{
	if (row != col) {
		return rg_Matrix();
	}

	rg_Matrix src( *this );

	rg_Matrix inverseMat(row, col);
	inverseMat.makeIdentity();

	rg_INT i, j, k;
	for (i=0; i<src.col; i++)  {
		k = -1;
		for ( j=i; j<src.row; j++)  {
			if ( rg_NZERO( src.element[j][i] ) )  {
				k = j;
				break;
			}
		}

		if ( k == -1 )  {
			//  no inverse matrix.
			return rg_Matrix();
		}

		//  i row <-> k row
		for (j=0; j<src.col; j++)  {
			rg_REAL temp1, temp2;
			temp1 = src.element[i][j];
			src.element[i][j] = src.element[k][j];
			src.element[k][j] = temp1;

			temp2 = inverseMat.element[i][j];
			inverseMat.element[i][j] = inverseMat.element[k][j];
			inverseMat.element[k][j] = temp2;
		}

		rg_REAL num = 1.0/src.element[i][i];
		for (j=0; j<src.col; j++)  {
			src.element[i][j] = src.element[i][j]*num;
			inverseMat.element[i][j] = inverseMat.element[i][j]*num;
		}
		for(k=0; k<row; k++)  {
			if (k!=i)
			{
				num = -1*src.element[k][i];
				for(j=0; j<src.col; j++)  {
					src.element[k][j] += num*src.element[i][j];
					inverseMat.element[k][j] += num*inverseMat.element[i][j];
				}
			}
		}
	}

	return inverseMat;
}

rg_Matrix rg_Matrix::getOneRowDeletedSubmatrix(const rg_INDEX& theCorrespondingRow) const
{
	rg_Matrix submatrix(row - 1, col);

	rg_INDEX index = 0;
	for(rg_INDEX i = 0;i < row;i++)
	{
		if(i != theCorrespondingRow)
		{
			for(rg_INDEX j = 0;j < col;j++)
			{
				submatrix[index][ j ] = element[ i ][ j ];
			}
			index++;
		}
	}
	
	return submatrix;
}

rg_Matrix rg_Matrix::getOneColDeletedSubmatrix(const rg_INDEX& theCorrespondingCol) const
{
	rg_Matrix submatrix(row, col - 1);

	rg_INDEX index = 0;
	for(rg_INDEX j = 0;j < col;j++)
	{
		if(j != theCorrespondingCol)
		{
			for(rg_INDEX i = 0;i < row;i++)
			{
				submatrix[ i ][index] = element[ i ][ j ];
			}
			index++;
		}
	}

	return submatrix;
}


rg_Matrix rg_Matrix::operator*(const rg_Matrix &matrix)const
{
    rg_Matrix mult(row, matrix.col);

    for(rg_INT i = 0; i < mult.row; i++)
        for(rg_INT j = 0; j < mult.col; j++ )
            for(rg_INT k = 0; k < col ; k++)
                mult.element[i][j] += element[i][k] * matrix.element[k][j]  ;

    return mult;
}

rg_Matrix rg_Matrix::operator*(const rg_REAL &n) const
{
    rg_Matrix output(row,col);
    for(rg_INT i = 0; i < row; i++)
        for(rg_INT j = 0; j < col; j++)
            output.element[i][j] = n*element[i][j];

    return output;
}

rg_Matrix rg_Matrix::operator/(const rg_REAL &n) const
{
    rg_Matrix dv(row,col);
    for(rg_INT i = 0; i < row; i++)
        for(rg_INT j = 0; j < col; j++)
            dv.element[i][j] = element[i][j]/n;

    return dv;
}
    
rg_Matrix rg_Matrix::operator+(const rg_Matrix &matrix) const
{
    rg_Matrix add;
    if ( col != matrix.col || row != matrix.row )
        return add;

    add.setSize(row,col);

    for(rg_INT i = 0; i < row; i++)
    {
        for(rg_INT j = 0; j < col; j++)
        {
            add.element[i][j] = element[i][j] + matrix.element[i][j];
        }
    }

    return add;
}

rg_Matrix rg_Matrix::operator-(const rg_Matrix &matrix) const
{
    rg_Matrix sub;
    if ( col != matrix.col || row != matrix.row )
        return sub;

    sub.setSize(row,col);

    for(rg_INT i = 0; i < row; i++)
        for(rg_INT j = 0; j < col; j++)
            sub.element[i][j] = element[i][j]-matrix.element[i][j];
    return sub;
}

rg_Matrix& rg_Matrix::operator=(const rg_Matrix &matrix)
{
    removeAll();
    setSize(matrix.row, matrix.col);

    for(rg_INT i = 0; i < row; i++)
        for(rg_INT j = 0; j < col; j++)
            element[i][j] = matrix.element[i][j];
    return *(this);
}

rg_REAL * rg_Matrix::operator[](const rg_INT& r)
{
    return element[r];
}

rg_Matrix operator*(const rg_REAL &scalar,const rg_Matrix& matrix)
{
    rg_Matrix mult(matrix.row, matrix.col);
    for(rg_INT i = 0; i < matrix.row; i++)
        for(rg_INT j = 0; j < matrix.col; j++)
            mult.element[i][j] = scalar*matrix.element[i][j];

    return mult;
}
/*
void rg_Matrix::show()
{
    for( rg_INT i=0 ;i<row; i++)
    {
        for( rg_INT j=0 ; j <col; j++)
        {
            cout.width(12);
            cout.precision(5);
            cout << element[i][j] << "  ";
        }
        cout << endl;
    }
}
*/
rg_REAL rg_Matrix::getNorm(const rg_INT &k) const
{
    rg_REAL maximum = 0.0;
    for(rg_INT i = 0; i < col; i++)
    {   
        rg_REAL sum = 0.0;
        for(rg_INT j = 0;j < row; j++)
            sum += pow(fabs(getElement(j, i)), k);
//          sum += pow(fabs(element[j][i]), k);

        sum = pow(sum, 1.0/k);  
        if( maximum < sum ) maximum = sum;
    }

    return maximum;
}

rg_REAL rg_Matrix::getNorm(const char *ch) const
{
    if( !strcmp(ch, "infinite"))
    {
        rg_Matrix A = transpose();
        return A.getNorm(1);
    }

    else 
    {
        return rg_NULL;
    }
}

// The following was inserted by Lee, Soon-Woong at 19 Jul. 1997.
// Two & Three-row determinants are inserted by Y. Cho at May 6, 2010.
rg_REAL rg_Matrix::determinant() const
{
    if(row != col) return rg_NULL;

    rg_REAL det = 1.0;
    if ( row == 2 ) {
        det = (element[0][0]*element[1][1]) - (element[0][1]*element[1][0]);
    }
    else if ( row == 3 ) {
        det =   (element[0][0]*element[1][1]*element[2][2]) 
              + (element[0][1]*element[1][2]*element[2][0]) 
              + (element[0][2]*element[1][0]*element[2][1]) 
              - (element[0][2]*element[1][1]*element[2][0]) 
              - (element[0][0]*element[1][2]*element[2][1]) 
              - (element[0][1]*element[1][0]*element[2][2]);
    }
    else {
        rg_Matrix m = *this;
        rg_Matrix b(row, col);
        fwdElimination(m, b);

        for(rg_INT i = 0; i < row; i++)
        {
            det *= m[i][i];
        }
    }

    return det;
}

void rg_Matrix::gaussianElimination(rg_Matrix &x, rg_Matrix &b)
{
    rg_Matrix m = *this;
    rg_INT rowOfb=b.getRowSize();
    rg_INT colOfb=b.getColSize();
    
    if (   x.getRowSize() != rowOfb 
        || x.getColSize() != colOfb  )
    {
        x.setSize( rowOfb,colOfb);
    }
    fwdElimination(m, b);
    bkSubstitution(m, x, b);
}

void rg_Matrix::fwdElimination(rg_Matrix &m, rg_Matrix &b) const
{
    rg_INT n = m.getRowSize();

    rg_REAL xmult = 0.0;

    for(rg_INT k = 0; k < n-1; k++)
    {
        for(rg_INT i = k+1; i < n; i++)
        {
/*            xmult = m[i][k]/m[k][k];
            m[i][k] = xmult;
            for(rg_INT j = k+1; j < n; j++)
            {
                m[i][j] = m[i][j] - xmult*m[k][j];
            }
            b[i][0] = b[i][0] - xmult*b[k][0];
			*/

            xmult = m[i][k]/m[k][k];
            m[i][k] = xmult;
            for(rg_INT j = k+1; j < n; j++)
            {
                m[i][j] = m[i][j] - xmult*m[k][j];
            }
            b[i][0] = b[i][0] - xmult*b[k][0];
        }
    }
}

void rg_Matrix::bkSubstitution(rg_Matrix &m, rg_Matrix &x, rg_Matrix &b) const
{
    if(row != col) return;
    
    rg_INT n = row;
    rg_REAL sum = 0.0;
    
    x[n-1][0]= b[n-1][0]/m[n-1][n-1];
    for(rg_INT i=n-2; i>=0; i--)
    {
        sum = b[i][0];
        for(rg_INT j=i+1; j<n; j++)
        {
            sum = sum - m[i][j]*x[j][0];
        }
        x[i][0] = sum/m[i][i];
    }
}

// LU decomposition carring out Doolittle's factorization
// Numerical Analysis pp. 131
// written by D. Kincaid and W. Cheney
void rg_Matrix::LUdecomposition(rg_Matrix &L, rg_Matrix &U)
{
    if(row != col) return;

//  rg_INT n = A.getRowSize();
    rg_INT n = row;
    for(rg_INT k=0; k<n; k++)
    {
        L[k][k] = 1.0;
        for(rg_INT j=k; j<n; j++)
        {
            rg_REAL sum = 0.0;
            for(rg_INT s=0; s<=k-1; s++)
                sum += L[k][s]*U[s][j];

            U[k][j] = element[k][j] - sum;
        }

        for(rg_INT i=k+1; i<n; i++)
        {
            rg_REAL sum = 0.0;
            for(rg_INT s=0; s<=k-1; s++)
                sum += L[i][s]*U[s][k];

            L[i][k] = (element[i][k] - sum)/U[k][k];
        }
    }
}
// End of Insertion

// 
rg_REAL*  rg_Matrix::getColMajorOrderedArray() const
{
    rg_REAL* output=rg_NULL;

    if ( element != rg_NULL )
    {        
        output=new rg_REAL[row*col];

        for( rg_INT j=0 ; j < col; j++ )
        {
            for( rg_INT i=0; i < row; i++)
            {
                output[ j*row+i ]=element[i][j];
            }
        }
    }

    return output;
}

rg_REAL*  rg_Matrix::getRowMajorOrderedArray() const
{

    rg_REAL* output=rg_NULL;

    if ( element != rg_NULL )
    {        
        output=new rg_REAL[row*col];

        for( rg_INT i=0 ; i < row; i++ )
        {
            for( rg_INT j=0; j < col; j++)
            {
                output[ i*col+j ]=element[i][j];
            }
        }
    }
  
    return output;
}

void rg_Matrix::display()
{

}
/*
std::ostream& operator<<(std::ostream& out, const rg_Matrix& mat)
{
	out << mat.row << "\t" << mat.col << endl;
	for(rg_INT i=0; i<mat.row; i++)
	{
		for(rg_INT j=0; j<mat.col; j++)
		{
			out << "" << mat.element[i][j] << '\t';
		}
		out << endl;
	}

	return out;

}

ofstream& operator<<(ofstream& fout, const rg_Matrix& mat)
{
	fout << mat.row << "\t" << mat.col << endl;
	for(rg_INT i=0; i<mat.row; i++)
	{
		for(rg_INT j=0; j<mat.col; j++)
		{
			fout << "" << mat.element[i][j] << '\t';
		}
		fout << endl;
	}

	return fout;

}
*/


