#ifndef _MANIFOLDIZEDBETASHAPE_H__
#define _MANIFOLDIZEDBETASHAPE_H__

#include "rg_Const.h"
#include "rg_dList.h"

#include "BetaUniverse.h"

using namespace V::GeometryTier;

#include "MBSBody.h"
#include "MBSShell.h"
#include "MBSFace.h"
#include "MBSEdge.h"
#include "MBSVertex.h"


//construction option
const rg_INT CONSTRUCT_WITH_WHOLE_BETA_SHAPE   = 0;
const rg_INT CONSTRUCT_WITH_REGULAR_BETA_SHAPE = 1;



class ManifoldizedBetaShape
{
private:
    BetaUniverse*           m_betaUniverse;
    rg_REAL                 m_betaValue;

    rg_dList< MBSBody* >    m_bodies;
    

public:
    ManifoldizedBetaShape();
    ~ManifoldizedBetaShape();

	
    //get function
    BetaUniverse*           getBetaUniverse() const;
    rg_REAL                 getBetaValue() const;
    rg_dList< MBSBody* >*   getBodies();
    
    rg_INT getNumOfBodies() const;
    rg_INT getNumOfShells() const;
    rg_INT getNumOfVertices() const;
    rg_INT getNumOfEdges() const;
    rg_INT getNumOfFaces()const;

    //  Youngsong Cho - 2011-04-04
    rg_INT getNumOfHandlesByIgnoringInteriorVoids() const;


    //set function
    void setBetaUniverse( BetaUniverse* betaUniverse );
    void setBetaValue( const rg_REAL& betaValue );
    void addBody( MBSBody* body );


    //  Youngsong Cho - 2011-04-04
    void construct( BetaUniverse* betaUniverse, const rg_REAL& betaValue, const rg_INT& constructOption = CONSTRUCT_WITH_WHOLE_BETA_SHAPE );



    void removeAll();
    rg_REAL getNumOfHoles(); //rg_REAL --> rg_INT 현재는 debugging을 위해
    
};

#endif
