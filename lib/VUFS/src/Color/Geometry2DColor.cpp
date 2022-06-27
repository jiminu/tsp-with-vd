#include "Geometry2DColor.h"

Geometry2DColor::Geometry2DColor()
{
    m_RGBForFaceColor     = DEFAULT_COLOR_OBJ[GRAY];
    m_RGBForBoundaryColor = DEFAULT_COLOR_OBJ[BLACK];
    m_RGBForCenterPTColor = DEFAULT_COLOR_OBJ[BLACK];
}

Geometry2DColor::Geometry2DColor(const Color4DrawingObj& RGBForFaceColor, const Color4DrawingObj& RGBForBoundaryColor, const Color4DrawingObj& RGBForCenterPTColor)
{
    m_RGBForFaceColor     = RGBForFaceColor;
    m_RGBForBoundaryColor = RGBForBoundaryColor;
    m_RGBForCenterPTColor = RGBForCenterPTColor;
}

Geometry2DColor::Geometry2DColor(const Geometry2DColor& geometry2DColor)
{
    copy_from(geometry2DColor);
}


Geometry2DColor& Geometry2DColor::operator=(const Geometry2DColor& geometry2DColor)
{
    if (this != &geometry2DColor)
    {
        copy_from(geometry2DColor);
    }

    return *this;
}


void Geometry2DColor::copy_from(const Geometry2DColor& geometry2DColor)
{
    m_RGBForFaceColor     = geometry2DColor.m_RGBForFaceColor;
    m_RGBForBoundaryColor = geometry2DColor.m_RGBForBoundaryColor;
    m_RGBForCenterPTColor = geometry2DColor.m_RGBForCenterPTColor;
}
