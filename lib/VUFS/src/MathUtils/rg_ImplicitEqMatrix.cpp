/////////////////////////////////////////////////////////////////////
//
//	  FILENAME    : rg_ImplicitEqMatrix.cpp
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

#include "rg_ImplicitEqMatrix.h"
#include <math.h>

rg_ImplicitEqMatrix::rg_ImplicitEqMatrix()
{
	col     = 0;
	row     = 0;
	element = rg_NULL;
}

/*
rg_ImplicitEqMatrix(const rg_DEGREE tdegree)
{
	rg_ImplicitEquation temp;
	for(rg_INT i = 0;i < tdegree;i++)
	{
		for(rg_INT j = 0;j < tdegree;j++)
		{

		}
	}
}
*/

rg_ImplicitEqMatrix::rg_ImplicitEqMatrix(const rg_INT &trow, const rg_INT &tcol)
{
	row = trow;
	col = tcol;

	element = new rg_ImplicitEquation*[row];

	for(rg_INT i = 0;i < row;i++)
	{
		element[ i ] = new rg_ImplicitEquation[col];
	}
}

rg_ImplicitEqMatrix::~rg_ImplicitEqMatrix()
{
	removeAll();				
}

rg_ImplicitEqMatrix::rg_ImplicitEqMatrix(const rg_ImplicitEqMatrix& SourceObj)
{
	row = SourceObj.row;
	col = SourceObj.col;
	
	element = new rg_ImplicitEquation*[row];
	rg_INT i = 0;
	for(i = 0;i < row;i++)
	{
		element[ i ] = new rg_ImplicitEquation[col];
	}

	for(i = 0;i < row;i++)
	{
		for(rg_INT j = 0;j < col;j++)
		{
			element[ i ][ j ] = SourceObj.getElement(i, j);
		}
	}
		
}

rg_INT rg_ImplicitEqMatrix::getRowSize() const
{
	return row;
}

rg_INT rg_ImplicitEqMatrix::getColSize() const
{
	return col;
}


rg_ImplicitEquation rg_ImplicitEqMatrix::getCofactor(const rg_INDEX& i, const rg_INDEX& j) const
{
	rg_ImplicitEqMatrix Rtemp = deleteOneRow( i );	
	rg_ImplicitEqMatrix Ctemp = Rtemp.deleteOneCol( j );

    rg_REAL scalar=pow(-1., i + j);
    rg_ImplicitEquation determinent=Ctemp.getDeterminant();
    rg_ImplicitEquation output=scalar*determinent;
	return  output ;
}


// implicitequaiton class operator*
// degree of output 

rg_ImplicitEquation rg_ImplicitEqMatrix::getDeterminant() const
{
	//show();//test
	//static rg_INT RowSize = row;
	//if(RowSize == 2)
	if(row == 2)
	{
		//rg_ImplicitEquation temp = element[ 0 ][ 0 ] * element[ 1 ][ 1 ]
		//	                   -element[ 0 ][ 1 ] * element[ 1 ][ 0 ];//test
		//temp.show();//test
		return element[ 0 ][ 0 ] * element[ 1 ][ 1 ]-element[ 0 ][ 1 ] * element[ 1 ][ 0 ];
	}
	else
	{
		// storage = 1(constant) + 2 dimensional * (degree) + Summation{i = degree}(i - 1)

		rg_ImplicitEquation output(1);
		rg_INT j = 0;

		//RowSize--;

		while(j < col)
		{			
			//rg_ImplicitEqMatrix Rtemp = deleteOneRow( 0 );//antique

			//Rtemp.show();//test
			
			//rg_ImplicitEqMatrix Ctemp = Rtemp.deleteOneCol( j );//antique

			rg_ImplicitEqMatrix temp = deleteOneElement(0, j);

			//Ctemp.show();//test
			//rg_ImplicitEquation temp = pow(-1, 0+1+j+1) * element[ 0 ][ j ] * Ctemp.getDeterminant();//test
			//temp.show();//test
			
			//output += (pow(-1, 0+1+j+1) * (element[ 0 ][ j ] * Ctemp.getDeterminant()));//antique
			//or output = output + pow(-1, 0+1+j+1) * element[ 0 ][ j ] * Ctemp.getDeterminant();//antique
            rg_REAL scalar=pow(-1., 0+1+j+1);
			output += ( scalar* (element[ 0 ][ j ] * temp.getDeterminant()));

			//output.show();//test
			j++;
		}
		return output;
	}	
}

rg_ImplicitEqMatrix rg_ImplicitEqMatrix::deleteOneRow(const rg_INDEX& theCorrespondingRow) const
{
	rg_ImplicitEqMatrix temp(row - 1, col);
	rg_INT i = 0, index = 0;

	while((i < row - 1) && (index < row))
	{
		if(index != theCorrespondingRow)
		{
			for(rg_INT j = 0;j < col;j++)
			{
				temp[ i ][ j ] = getElement(index, j);
			}
			i++;
		}
		index++;
	}

	return temp;
}

rg_ImplicitEqMatrix rg_ImplicitEqMatrix::deleteOneCol(const rg_INDEX& theCorrespondingCol) const
{
	rg_ImplicitEqMatrix temp(row, col - 1);
	rg_INT j = 0, index = 0;

	while((j < col - 1) && (index < col))
	{
		if(index != theCorrespondingCol)
		{
			for(rg_INT i = 0;i < row;i++)
			{
				temp[ i ][ j ] = getElement(i, index);
			}
			j++;
		}
		index++;		
	}

	return temp;	
}

rg_ImplicitEqMatrix rg_ImplicitEqMatrix::deleteOneElement(const rg_INDEX& theCorrespondingRow, const rg_INDEX& theCorrespondingCol) const
{
	rg_ImplicitEqMatrix oneRowDeleted(row - 1, col);
	rg_INT i = 0, j = 0, index = 0;

	while((i < row - 1) && (index < row))
	{
		if(index != theCorrespondingRow)
		{
			for(j = 0;j < col;j++)
			{
				oneRowDeleted[ i ][ j ] = getElement(index, j);
			}
			i++;
		}
		index++;
	}

	j = index = 0;

	rg_ImplicitEqMatrix oneElementDeleted(row - 1, col - 1);

	while((j < col - 1) && (index < col))
	{
		if(index != theCorrespondingCol)
		{
			for(i = 0;i < row - 1;i++)
			{
				oneElementDeleted[ i ][ j ] = oneRowDeleted.getElement(i, index);
			}
			j++;
		}
		index++;		
	}

	return oneElementDeleted;
}

void rg_ImplicitEqMatrix::removeAll()
{
	if(element != rg_NULL)
	{
		for(rg_INT i = 0;i < row;i++)
			delete[] element[ i ];
		delete []element;
	}
		
	element = rg_NULL;
	row = 0;
	col = 0;
}

rg_ImplicitEquation rg_ImplicitEqMatrix::getElement(const rg_INDEX& i, const rg_INDEX& j) const
{
	return element[ i ][ j ];
}

rg_ImplicitEquation* rg_ImplicitEqMatrix::operator[](const rg_INDEX& index)
{
	return element[ index ];
}

rg_ImplicitEqMatrix& rg_ImplicitEqMatrix::operator=(const rg_ImplicitEqMatrix& SourceObj)
{
	if(this == &SourceObj)
		return *this;

	removeAll();

	row = SourceObj.row;
	col = SourceObj.col;
	
	element = new rg_ImplicitEquation*[row];
	rg_INT i = 0;
	for(i = 0;i < row;i++)
	{
		element[ i ] = new rg_ImplicitEquation[col];
	}

	for(i = 0;i < row;i++)
	{
		for(rg_INT j = 0;j < col;j++)
		{
			element[ i ][ j ] = SourceObj.getElement(i, j);
		}
	}

	return *this;
}

void rg_ImplicitEqMatrix::show() const
{
	/*
	for(rg_INT i = 0;i < row;i++)
	{
		for(rg_INT j = 0;j < col;j++)
		{
			cout << i << j << '\n';
			element[ i ][ j ].show();
		}
	}
	*/
}



