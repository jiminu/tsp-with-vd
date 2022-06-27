#include <math.h>

#include "rg_Const.h"
//#include "..\header\rg_Const.h"
//#include "..\header\pt2D.h"
//#include "..\header\vector2D.h"

// each of the following functions use a default argument: res=rg_MATH_RES 
// since the default argument should NOT be redefined
// the #ifndef ... #define ... #endif
// have been used in the ???h file.

rg_REAL rg_ABS(const rg_REAL& x)
{
   return (x>0 ? x : -x);
}

// relative operators
rg_INT rg_EQ(const rg_REAL& r1, const rg_REAL& r2, const rg_REAL& res)
        { return (rg_ABS(r1 - r2) <= res); }
rg_INT rg_GT(const rg_REAL& r1, const rg_REAL& r2, const rg_REAL& res)
        { return (r1 > r2 + res); }
rg_INT rg_LT(const rg_REAL& r1, const rg_REAL& r2, const rg_REAL& res)
        { return (r1 < r2 - res); }
rg_INT rg_GE(const rg_REAL& r1, const rg_REAL& r2, const rg_REAL& res)
        { return (! rg_LT(r1, r2, res)); }
rg_INT rg_LE(const rg_REAL& r1, const rg_REAL& r2, const rg_REAL& res)
        { return (! rg_GT(r1, r2, res)); }
rg_INT rg_NE(const rg_REAL& r1, const rg_REAL& r2, const rg_REAL& res)
        { return (! rg_EQ(r1, r2, res)); }

//      Compare to rg_ZERO        
rg_INT rg_ZERO(const rg_REAL& r, const rg_REAL& res) 
		{ return (rg_ABS(r) <= res); }
rg_INT rg_POS(const rg_REAL& r, const rg_REAL& res) 
		{ return (r > res); }
rg_INT rg_NEG(const rg_REAL& r, const rg_REAL& res) 
		{ return (r < -res); }
rg_INT rg_NZERO(const rg_REAL& r, const rg_REAL& res) 
		{ return (! rg_ZERO(r, res)); }
rg_INT rg_NPOS(const rg_REAL& r, const rg_REAL& res) 
		{ return (! rg_POS(r, res)); }
rg_INT rg_NNEG(const rg_REAL& r, const rg_REAL& res) 
		{ return (! rg_NEG(r, res)); }

//      Between Two Ordered Reals: (the real vaues are ordered, i.e., left <= right)    
rg_INT rg_BTOR(const rg_REAL& left, const rg_REAL& value, const rg_REAL& right, const rg_REAL& res)
		{ return (rg_LE(left, value, res) && rg_LE(value, right, res)); }
rg_INT rg_BTORexclusive(const rg_REAL& left, const rg_REAL& value, const rg_REAL& right, const rg_REAL& res)
		{ return (rg_LT(left, value, res) && rg_LT(value, right, res)); }

//      Between Two Reals: (the real values are not ordered.)    
rg_INT rg_BTR(const rg_REAL& r1, const rg_REAL& r, const rg_REAL& r2, const rg_REAL& res)
		{ return ((rg_LE(r1, r, res) && rg_LE(r, r2, res)) || (rg_LE(r2, r, res) && rg_LE(r, r1, res))); }
 
//----- if r1 > r2 then 1 is returned  
//----- if r1 = r2 then 0 is returned
//----- if r1 < r2 then -1is  returned 
rg_INT rg_COMPARE(const rg_REAL& r1, const rg_REAL& r2, const rg_REAL& res)
{
   if ( rg_ABS(r1 - r2) < res ) 
      return 0;
   else if ( r1 > r2 )      
      return 1;
   else 
      return -1;
}


rg_REAL   rg_Max( rg_REAL const& first, rg_REAL const& second)
{ if( first > second ) return first;
	else return second;
}
rg_REAL   rg_Min( rg_REAL const& first, rg_REAL const& second)
{ if( first < second ) return first;
	else return second;
}


