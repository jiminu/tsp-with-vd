#ifndef _VIDICITEM_H
#define _VIDICITEM_H

#include "ConstForVoronoiDiagram3D.h"

#include "rg_Point3D.h"

#include "rg_BallGenerator.h"
#include "VDVertex.h"

namespace V {

namespace GeometryTier {


class VIDIC;

class VIDICItem // Vertex Index DICtionary
{
private:
    BallGenerator*   m_balls[NUM_DEFINING_CELLS_OF_ONE_VERTEX];
    VDVertex*        m_vertex;
    rg_FLAG          m_AvailableNumVertexByBallConfiguration;

public:

    friend class VIDIC;

    //  constructor & deconstructor..
    VIDICItem();
    VIDICItem( BallGenerator** balls, VDVertex* vertex );
    VIDICItem( BallGenerator** balls, VDVertex* vertex, const rg_FLAG& numVertexBy4Balls );
    VIDICItem( BallGenerator* ball1, BallGenerator* ball2, BallGenerator* ball3, BallGenerator* ball4,VDVertex* vertex );
    VIDICItem( const VIDICItem& anItem );
    ~VIDICItem();

    //  get functions..     
    BallGenerator*   getBall(const rg_INT& i) const;
    BallGenerator**  getAllBalls();

    VDVertex* getVertex() const;

    rg_FLAG  isSameVertexConfiguration( BallGenerator* ball1, 
                                        BallGenerator* ball2, 
                                        BallGenerator* ball3, 
                                        BallGenerator* ball4, 
                                        const rg_Point3D& vertexCoord );

    //  set functions..
    void    setBall( const rg_INT& i, BallGenerator* aBall );
    void    setAllBalls( BallGenerator** balls );
    rg_FLAG setBall( BallGenerator* aBall );

    void setVertex( VDVertex* aVertex );

    //  operator overloading..
    VIDICItem& operator =( const VIDICItem& anItem );
};

} // namespace GeometryTier

} // namespace V

#endif
