#ifndef _RG_RANDOMNUMBER_H
#define _RG_RANDOMNUMBER_H

#include <time.h>
#include "rg_Const.h"
const unsigned long maxshort   = 65536L;
const unsigned long multiplier = 1194211693L;
const unsigned long adder      = 12345L;

class rg_RandomNumber
{
private:
	unsigned long randSeed;

public:
	rg_RandomNumber();
	rg_RandomNumber(unsigned long s);
	~rg_RandomNumber();

	void setRandomSeed(unsigned long s);
	unsigned short random(unsigned long n);
	rg_REAL fRandom();
};

#endif


