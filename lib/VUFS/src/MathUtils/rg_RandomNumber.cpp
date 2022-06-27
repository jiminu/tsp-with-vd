#include "rg_RandomNumber.h"

rg_RandomNumber::rg_RandomNumber()
{
}

rg_RandomNumber::rg_RandomNumber(unsigned long s = 0)
{
	if(s == 0)
		randSeed = time(0);
	else
		randSeed = s;
}

rg_RandomNumber::~rg_RandomNumber()
{
}

void rg_RandomNumber::setRandomSeed(unsigned long s)
{
	randSeed = s;
}

unsigned short rg_RandomNumber::random(unsigned long n)
{
	randSeed = multiplier * randSeed + adder;
	return (unsigned short)((randSeed >> 16) % n);
}

rg_REAL rg_RandomNumber::fRandom()
{
	return random(maxshort)/rg_REAL(maxshort);
}	


