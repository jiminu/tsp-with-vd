#include "FunctionsForPocket.h"

#include "Pocket.h"
#include "PocketPrimitive.h"



#include "rg_Const.h"
#include "rg_RelativeOp.h"


int comparePocketPrimitivesByNumFaces(const void* pocket1, const void* pocket2)
{
    PocketPrimitive* primitive1 = *(PocketPrimitive**) pocket1;
    PocketPrimitive* primitive2 = *(PocketPrimitive**) pocket2;

    rg_INT numFacesOfPrimitive1 = primitive1->getNumMBSFaces();
    rg_INT numFacesOfPrimitive2 = primitive2->getNumMBSFaces();

    if ( numFacesOfPrimitive1 > numFacesOfPrimitive2 )
        return -1;
    else if ( numFacesOfPrimitive1 == numFacesOfPrimitive2 )
        return 0;
    else 
        return 1;
}



int comparePocketPrimitivesByAreaOfFaces(const void* pocket1, const void* pocket2)
{
    PocketPrimitive* primitive1 = *(PocketPrimitive**) pocket1;
    PocketPrimitive* primitive2 = *(PocketPrimitive**) pocket2;

    rg_REAL areaOfFacesOfPrimitive1 = primitive1->getAreaOfMBSFaces();
    rg_REAL areaOfFacesOfPrimitive2 = primitive2->getAreaOfMBSFaces();

    if ( areaOfFacesOfPrimitive1 > areaOfFacesOfPrimitive2 )
        return -1;
    else if ( areaOfFacesOfPrimitive1 == areaOfFacesOfPrimitive2 )
        return 0;
    else 
        return 1;
}


