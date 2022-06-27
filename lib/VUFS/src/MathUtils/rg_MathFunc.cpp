#include "rg_MathFunc.h"
#include "rg_RelativeOp.h"

#include <time.h>
#include <stdlib.h>

rg_REAL rg_MathFunc::cubicRoot(const rg_REAL &r)
{
	if(rg_NNEG(r, rg_SYSTEM_RES))
		return pow(r, 1.0/3.0);
	else
		return -pow(-r, 1.0/3.0);
}

rg_ComplexNumber *rg_MathFunc::cubicRoot(const rg_ComplexNumber& c)
{
	return rg_MathFunc::root(c, 3);
}

rg_ComplexNumber *rg_MathFunc::root(const rg_ComplexNumber& c, const rg_INT& n)
{
	rg_REAL a = c.getRealNumber();
	rg_REAL b = c.getImaginaryNumber();
	rg_REAL r = sqrt( pow(a, 2) + pow(b, 2));
	rg_REAL PIE = 4.0 * atan(1.0);

	// need to handle cases that a == 0 && b == 0.
	if( c.isZero() ) 
	{
		rg_ComplexNumber *root = new rg_ComplexNumber[n];
		for(rg_INT i=0; i< n; i++)
		{
			root[i] = rg_ComplexNumber(0, 0);
		}
		return root;
	}

   	rg_INT  flag = 0;
	if( rg_POS(a)  && rg_NNEG(b) ) flag = 1;
	if( rg_NEG(a)  && rg_NNEG(b) ) flag = 2;
	if( rg_NEG(a)  && rg_NEG(b)  ) flag = 3;
	if( rg_NNEG(a) && rg_NEG(b)  ) flag = 4;

	rg_REAL theta = 0;
	switch(flag)
	{
		case 1: theta = acos(a/r); break;
		case 2: theta = acos(a/r); break;
		case 3: theta = 2 * PIE - acos(a/r); break;
		case 4: theta = 2 * PIE - acos(a/r); break;
	}


	rg_REAL rhou = pow(r, 1.0/(rg_REAL)n);
	
	rg_ComplexNumber *root = new rg_ComplexNumber[n];
	for(rg_INT i=0; i< n; i++)
	{
		root[i] = rg_ComplexNumber(rhou * cos(theta/n + 2*i*PIE/n), 
			                    rhou * sin(theta/n + 2*i*PIE/n) );
	}

/*
// We test the value resulted from above.
	rg_INT okCount = 0;
	for(i = 0; i < n; i++)
	{
		rg_ComplexNumber vTest;
		for(rg_INT j=0; j < n; j++)
		{
			vTest = vTest * root[i];
		}
		if(vTest.isZero()) okCount++;
	}

	if(okCount == n)
		cout << "We found " << n << "th root correctly!!" << endl;
*/
	return root;
}

/*
rg_Point2D *rg_MathFunc::solveQuadraticEq(const rg_Point2D &a, const rg_Point2D &b, const rg_Point2D &c)
{
	rg_Point2D t = rg_Point2D(sqrt(b.getX()*b.getX() - 4*a.getX()*c.getX()),
		                sqrt(b.getY()*b.getY() - 4*a.getY()*c.getY()) );
	
	rg_Point2D *solution = new rg_Point2D[2]; 
	for(rg_INT i=0; i < 2; i++)
	{
		solution[i].setX((-b.getX() + pow(-1,i)*t.getX())/2.0/a.getX());
		solution[i].setY((-b.getY() + pow(-1,i)*t.getY())/2.0/a.getY());
	}

	return solution;
}
*/

// solution of an equation, a x^2 + b x + c
rg_REAL *rg_MathFunc::solveQuadraticEq(const rg_REAL &a, const rg_REAL &b, const rg_REAL &c)
{
	rg_REAL *solution = new rg_REAL[2]; 

	// a의 값이 Zero인 경우를 생각해야 한다. 
	if( rg_ZERO(a) ) 
	{
		solution[0] = -c/b;
		solution[1] = rg_NULL;
	}

	else
	{
		for(rg_INT i=0; i < 2; i++)
			solution[i] = (-b + pow(-1.,i)*sqrt(b*b-4.*a*c))/(2.0*a);
	}
	return solution;
}

rg_ComplexNumber *rg_MathFunc::solveQuadraticEq(rg_ComplexNumber &a, rg_ComplexNumber &b, rg_ComplexNumber &c, rg_INT& num)
{
	rg_ComplexNumber *solution = rg_NULL;
	if(a.isPureRealNumber() && rg_ZERO(a.getRealNumber()))
	{
		solution = new rg_ComplexNumber[ 1 ];
		solution[ 0 ] = -c / b;
		num = 1;

	}
	else
	{
		solution = new rg_ComplexNumber[ 2 ];
		rg_ComplexNumber* sRoot = root(b*b-4*a*c, 2);

		rg_ComplexNumber tSol[ 4 ];
		
		tSol[ 0 ] = (-b + sRoot[ 0 ])/(2.0*a);
		tSol[ 1 ] = (-b - sRoot[ 0 ])/(2.0*a);
		tSol[ 2 ] = (-b + sRoot[ 1 ])/(2.0*a);
		tSol[ 3 ] = (-b - sRoot[ 1 ])/(2.0*a);

		solution[ 0 ] = tSol[ 0 ];
		for(rg_INT i = 1;i <= 3;i++)
		{
			if(!(solution[ 0 ] == tSol[ i ]))
			{
				solution[ 1 ] = tSol[ i ];
				break;
			}
		}
		num = 2;
		delete [] sRoot;
	}
	return solution;
}

rg_REAL rg_MathFunc::minValue(const rg_REAL &a, const rg_REAL &b, const rg_REAL &c)
{
	rg_REAL tm_min = a;

	if( rg_GT(tm_min, b) )
	{
		tm_min = b;
	}
	else if( rg_GT(tm_min, c) )
	{
		tm_min = c;
	}

	return tm_min;
}

rg_REAL rg_MathFunc::maxValue(const rg_REAL &a, const rg_REAL &b, const rg_REAL &c)
{
	rg_REAL tm_max = a;

	if( rg_LT(tm_max, b) )
	{
		tm_max = b;
	}
	else if( rg_LT(tm_max, c) )
	{
		tm_max = c;
	}
	
	return tm_max;
}
//
//  
//        n*(n-1)* ... *(n-i+1)
// nCi = ---------------------------
//        i*(i-1)* ... *1
rg_INT  rg_MathFunc::combination(const rg_INT& n, const rg_INT& i)
{
/*
	if( i < 0 || i > n || n < 0 )
	{	
		return 0.0;
	}
	if ( i == 0 || i == n )
	{
		return 1.0;
	}
*/

    rg_INT x=( i > n/2 ) ? n-i: i;// nCi = nCn-i
                               // and costs small computation in min(i,n-i) case

	rg_REAL output=1.0;
	for( rg_INT j=0; j < x ; j++ )
	{
		output*=( (rg_REAL)(n-j) )/(x-j);
	}

	return (rg_INT)output;
}

rg_REAL rg_MathFunc::factorial(const rg_INT &n)
{
	rg_REAL output=1.0;
	
	for( rg_INT i=2 ; i <= n; i++ )
	{
		output*=i;
	}

	return output;
}



void rg_MathFunc::compute_mean_and_standard_deviation(const list<rg_REAL>& observations, rg_REAL & mean, rg_REAL & standardDeviation)
{
    //rg_REAL sum = accumulate(observations.begin(), observations.end(), 0.0);
    rg_REAL summation = 0.0;
    list<rg_REAL>::const_iterator i_observation = observations.begin();
    for (; i_observation != observations.end(); ++i_observation)
    {
        summation += (*i_observation);
    }

    mean = summation / observations.size();

    rg_REAL summationOfSquaredDeviations = 0.0;
    i_observation = observations.begin();
    for (; i_observation != observations.end(); ++i_observation)
    {
        rg_REAL deviation = (*i_observation) - mean;
        summationOfSquaredDeviations += (deviation * deviation);
    }

    standardDeviation = sqrt(summationOfSquaredDeviations / observations.size());
}



// Generating combinatorial objects
// Generating all the possible subset of 'sizeOfSet'
// Input set must be of the form '0, 1, 2, ...'
rg_INT** rg_MathFunc::enumerateSubset(rg_INT* set, rg_INT sizeOfSet, rg_INT sizeOfSubset)
{
	rg_INT noOfSubset = rg_MathFunc::combination(sizeOfSet, sizeOfSubset);

	rg_INT** allPossibleSubsets = new rg_INT* [noOfSubset];
	rg_INT i = 0;
	for(i = 0;i < noOfSubset;i++)
	{
		allPossibleSubsets[ i ] = new rg_INT [sizeOfSubset];
	}

	// enumerating k-combinations from n-set

	// Initialization
	rg_INT indexOfSubsets = 0;
	rg_INT* buffer = new rg_INT[sizeOfSubset];
	for(i = 0;i < sizeOfSubset;i++)
	{
		allPossibleSubsets[indexOfSubsets][ i ] = buffer[ i ] = set[ i ];
	}

	// Generating all subsets
	indexOfSubsets++;
	while(indexOfSubsets < noOfSubset)
	{
		rg_INT indexOfBuffer = 0, indexOflargest = 0;;
		while(indexOfBuffer < sizeOfSubset)
		{
			allPossibleSubsets[indexOfSubsets][indexOfBuffer] = buffer[indexOfBuffer];

			if(buffer[indexOfBuffer] < sizeOfSet - sizeOfSubset + indexOfBuffer)
			{
				indexOflargest = indexOfBuffer;
			}
			indexOfBuffer++;
		}

		indexOfBuffer = indexOflargest;
		rg_INT largest = buffer[indexOflargest];

		while(indexOfBuffer < sizeOfSubset)
		{
			allPossibleSubsets[indexOfSubsets][indexOfBuffer] = buffer[indexOfBuffer] = largest + 1;
			largest++;
			indexOfBuffer++;
		}

		indexOfSubsets++;
	}
	delete[] buffer;
	return allPossibleSubsets;
}

// Generating combinatorial objects
// Generating all the possible sequences that consists of 0's and 1's
rg_INT** rg_MathFunc::enumerateZeroOneSequence(rg_INT noOfZero, rg_INT noOfOne)
{
	// initialization
	rg_INT sizeOfSet = noOfZero + noOfOne;
	rg_INT* set = new rg_INT [sizeOfSet];
	rg_INT sizeOfSubset = noOfOne;
	//rg_INT* subset = new rg_INT [sizeOfSubset];

	rg_INT i = 0;
	for(i = 0;i < sizeOfSet;i++)
		set[ i ] = i;

	// enumerating all the possible positions of one's 
	//             and all the possible zero-one sequences
	rg_INT noOfPossibleSubset = rg_MathFunc::combination(sizeOfSet, sizeOfSubset);

	// Initialization
	rg_INT** allPossibleZeroOneSequences = new rg_INT* [noOfPossibleSubset];
	for(i = 0;i < noOfPossibleSubset;i++)
	{
		allPossibleZeroOneSequences[ i ] = new rg_INT[sizeOfSet];
		for(rg_INT j = 0;j < sizeOfSet;j++)
			allPossibleZeroOneSequences[ i ][ j ] = 0;
	}

	rg_INT indexOfSubsets = 0;
	rg_INT* buffer = new rg_INT[sizeOfSubset];
	for(i = 0;i < sizeOfSubset;i++)
	{
		buffer[ i ] = set[ i ];
		allPossibleZeroOneSequences[indexOfSubsets][buffer[ i ]] = 1;
	}
	delete[] set;
	// Generating all subsets
	indexOfSubsets = 1;
	while(indexOfSubsets < noOfPossibleSubset)
	{
		rg_INT indexOfBuffer = 0, indexOflargest = 0;;
		while(indexOfBuffer < sizeOfSubset)
		{
			if(buffer[indexOfBuffer] < sizeOfSet - sizeOfSubset + indexOfBuffer)
			{
				indexOflargest = indexOfBuffer;
			}
			indexOfBuffer++;
		}
		indexOfBuffer = indexOflargest;
		rg_INT largest = buffer[indexOflargest];
		while(indexOfBuffer < sizeOfSubset)
		{
			buffer[indexOfBuffer] = largest + 1;
			largest++;
			indexOfBuffer++;
		}

		rg_INT i = 0;
		while(i < sizeOfSubset)
		{
			allPossibleZeroOneSequences[indexOfSubsets][buffer[ i ]] = 1;
			i++;
		}
		indexOfSubsets++;
	}
	return allPossibleZeroOneSequences;
}

// Generating combinatorial objects
// Generating all the possible sequences that consists of 0's and 1's
rg_INT** rg_MathFunc::enumerateZeroOneSequenceRevised(rg_INT noOfZero, rg_INT noOfOne)
{
	// Initialization
	rg_INT sizeOfOnePath = noOfZero + noOfOne;

	rg_INT noOfAllPaths = rg_MathFunc::combination(sizeOfOnePath, noOfOne);

	rg_INT** allPossibleZeroOneSequences = new rg_INT* [noOfAllPaths];
    rg_INT** index = new rg_INT* [noOfAllPaths];
	rg_INT i = 0;
	for(i = 0;i < noOfAllPaths;i++)
	{
		allPossibleZeroOneSequences[ i ] = new rg_INT[sizeOfOnePath];
		index[ i ] = new rg_INT[noOfOne];
		rg_INT j = 0;
		for(j = 0;j < sizeOfOnePath;j++)
			allPossibleZeroOneSequences[ i ][ j ] = 0;
		if(i == 0)
		{
			for(j = 0;j < noOfOne;j++)
				index[ i ][ j ] = j + 1;
		}
	}

	// enumerating all the possible positions of one's 
	for(i = 1;i < noOfAllPaths;i++)
	{
		for(rg_INT j = 0;j < noOfOne;j++)
		{
			if(index[i - 1][noOfOne - j - 1] < sizeOfOnePath - j)
			{
				rg_INT k = 0;
				for(k = 0;k < noOfOne - j - 1;k++)
					index[ i ][ k ] = index[i - 1][ k ];
				index[ i ][noOfOne - j - 1] = index[i - 1][noOfOne - j - 1] + 1;

				for(k = noOfOne - j;k < noOfOne;k++)
					index[ i ][ k ] = index[ i ][k - 1] + 1;
				break;
			}
		}
	}

	// enumerating all the possible zero-one sequences(to be deleted)
	for(i = 0;i < noOfAllPaths;i++)
	{
		for(rg_INT j = 0;j < noOfOne;j++)
		{
			allPossibleZeroOneSequences[ i ][index[ i ][ j ] - 1] = 1;
		}
	}
	return allPossibleZeroOneSequences;
}


rg_REAL rg_MathFunc::randomNumber(const rg_REAL& min, const rg_REAL& max)
{
    // rand() function returns a pseudorandom integer in the range 0 to RAND_MAX.
    // we should use srand() function to seed the pseudorandom-number generator before calling rand.

    return min + (double)rand() / RAND_MAX * (max - min);
}

void   rg_MathFunc::shuffleSequence(rg_INT*  inputSequence,
                                    rg_INT*& shuffledSequence, 
                                    const rg_INT& seqSize)
{
    // initialize boolean array which is used for storing whether each number has occurred in the previous run.
    rg_INDEX i;
    rg_BOOL* bAlreadyOccurred = rg_NULL;
    bAlreadyOccurred = new rg_BOOL[seqSize];
    for (i = 0;i < seqSize;i++)
    {
        bAlreadyOccurred[ i ] = rg_FALSE;
    }

    srand((unsigned)time(NULL));

    // generate random number in a sequence which does not occur in the previous run
    rg_INDEX j = 0;
    for (i = seqSize - 1;i >= 0;i--)
    {
        rg_INDEX currRandomIndex = randomNumber(0, i);

        rg_INDEX k = 0;
        rg_INDEX currIndex = 0;
        while (k <= currRandomIndex)
        {
            if(! bAlreadyOccurred[ currIndex++ ])
                k++;
        }
        
        shuffledSequence[ j++ ] = inputSequence[currIndex - 1];
        bAlreadyOccurred[currIndex - 1] = rg_TRUE;
    }

    if(bAlreadyOccurred != rg_NULL)
        delete [] bAlreadyOccurred;
}
rg_REAL rg_MathFunc::compute_root_of_polynomial_via_Newton_method_with_initial_solution(rg_Polynomial & polynomial, const rg_REAL& initialSolution, const rg_REAL& terminatingTolerance /*= rg_MATH_RES*/, const rg_INT& numIterations /*= 40*/)
{
	rg_REAL prevSol = initialSolution;
	rg_REAL currSol = initialSolution;
	rg_REAL funcValue = polynomial.evaluatePolynomial(prevSol);
	rg_Polynomial derivativeOfPolynomial = polynomial.makeDerivative();
	
	rg_INT i;
	for (i = 0; i < numIterations; ++i)
	{
		rg_REAL derivative = derivativeOfPolynomial.evaluatePolynomial(prevSol);
		currSol = prevSol - funcValue / derivative;
		funcValue = polynomial.evaluatePolynomial(currSol);
		
		if(rg_EQ(rg_ABS(currSol-prevSol), terminatingTolerance) || rg_ZERO(funcValue, terminatingTolerance))
			break;

		prevSol = currSol;
	}

	//std::cout << "initial solution: " << initialSolution << std::endl;
	//std::cout << "#iteration: " << i << std::endl;
	//std::cout << "solution: " << currSol << std::endl;

	return currSol;
}
