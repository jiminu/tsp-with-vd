#ifndef _GENERATOR_H
#define _GENERATOR_H

#include "rg_Const.h"

#include "Sphere.h"

namespace V {

namespace GeometryTier {



const rg_INT  NO_TYPE   = 0;
const rg_INT  BALL_TYPE = 1;

class VDCell;

class Generator
{
protected:
    rg_INT  m_ID;
    VDCell* m_cell;

    rg_FLAG m_ProbeTangibility; 


public:
    rg_FLAG m_checkForBucket;

public:
    //  constructor & deconstructor..
    Generator();
    Generator(const rg_INT& ID);
    Generator(const rg_INT& ID, VDCell* cell);
    Generator(const Generator& aGenerator);
    virtual ~Generator();

    //  get functions.. 
    rg_INT  getID() const;
    VDCell* getCell() const;

    virtual rg_INT  getTypeOfGenerator() const;
    virtual rg_FLAG isThisGeneratorIntersectedWith(Generator* anotherGenerator);
    virtual rg_FLAG isThisGeneratorIntersectedWith(const Sphere& aSphere);

    virtual Generator* duplicate() const;

    rg_FLAG    isTangible() const;
    void       isTangible(const rg_FLAG& tangibility);

    //  set functions..
    void setID(const rg_INT& ID);
    void setCell(VDCell* cell);
    void setGenerator(const rg_INT& ID, VDCell* cell);
    void setGenerator(const Generator& aGenerator);

    //  operator overloading..
    Generator& operator =(const Generator& aGenerator);

};

} // namespace GeometryTier

} // namespace V


#endif


