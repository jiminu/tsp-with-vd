#ifndef _BALL_H
#define _BALL_H

#include "Sphere.h"
#include "TopologicalEntity.h"

namespace V {

namespace GeometryTier {


class Ball : public TopologicalEntity
{
private: 
    Sphere m_geometry;
    void*  m_property;
    
    rg_INT m_IDFromInput;

public:
    Ball();
    Ball(const rg_INT& ID);
    Ball(const Sphere& sphere, 
         void* property = rg_NULL, 
         const rg_INT& IDFromInput = -1);
    Ball(const rg_INT& ID, const Sphere& sphere, 
         void* property = rg_NULL, 
         const rg_INT& IDFromInput = -1);
    Ball(const Ball& ball);
    ~Ball();

    Sphere getGeometry() const;
    void*  getProperty() const;
    rg_INT getIDFromInput() const;

    void   setGeometry(const Sphere& sphere);
    void   setProperty(void* property);
    void   setIDFromInput(const rg_INT& IDFromInput);

    Ball&  operator=(const Ball& ball);
};

} // namespace GeometryTier
} // namespace V


#endif 


