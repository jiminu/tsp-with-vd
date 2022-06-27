#ifndef _BALLGENERATOR_H
#define _BALLGENERATOR_H

#include "rg_Const.h"

#include "rg_Point3D.h"

//#include "Generator.h"
#include "Sphere.h"
//#include "Attribute.h"

namespace V {

namespace GeometryTier {


class VDCell;

class BallGenerator : /*public Generator ,*/ public Sphere
{
private:
    rg_INT   m_ID;
    VDCell*  m_cell;

    void*    m_property;
    rg_INT   m_IDFromInput;


    rg_FLAG  m_ProbeTangibility; 

public:
    rg_FLAG m_checkForBucket;


    //  constructor & deconstructor..
    BallGenerator();
    BallGenerator(const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, const rg_REAL& radius, void* property = rg_NULL, const rg_INT& IDFromInput = -1);
    BallGenerator(const rg_Point3D& center, const rg_REAL& radius, void*  property = rg_NULL, const rg_INT& IDFromInput = -1);
    BallGenerator(const rg_INT& ID, const rg_REAL& x, const rg_REAL& y, const rg_REAL& z, const rg_REAL& radius, void* property = rg_NULL, const rg_INT& IDFromInput = -1);
    BallGenerator(const rg_INT& ID, const rg_Point3D& center, const rg_REAL& radius, void* property = rg_NULL, const rg_INT& IDFromInput = -1);
    BallGenerator(const rg_INT& ID,  const Sphere& ball, void* property = rg_NULL, const rg_INT& IDFromInput = -1);
    BallGenerator(const BallGenerator& aBallGenerator);
    ~BallGenerator();

    //  get functions.. 
    rg_INT     getID() const;
    void*      getProperty() const;
    VDCell*    getCell() const;
    Sphere     getBall() const;
    rg_INT     getIDFromInput() const;

    rg_FLAG    isThisGeneratorIntersectedWith(BallGenerator* anotherGenerator);
    rg_FLAG    isThisGeneratorIntersectedWith(const Sphere& aSphere);

    //Generator* duplicate() const;

    rg_FLAG    isTangible() const;
    void       isTangible(const rg_FLAG& tangibility);

    //  set functions..
    void setID(const rg_INT& ID);
    void setCell(VDCell* cell);
    void setProperty(void* property);
    void setIDFromInput(const rg_INT& IDFromInput);

    //  operator overloading..
    BallGenerator& operator =(const BallGenerator& aBallGenerator);
};

} // namespace GeometryTier

} // namespace V


#endif


