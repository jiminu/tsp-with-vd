#ifndef _2D_GEOMETRY_COLOR_
#define _2D_GEOMETRY_COLOR_

#include "ConstForColor.h"

class Geometry2DColor
{
private:
    Color4DrawingObj  m_RGBForFaceColor;
    Color4DrawingObj  m_RGBForBoundaryColor;
    Color4DrawingObj  m_RGBForCenterPTColor;

public:
    Geometry2DColor();
    Geometry2DColor(const Color4DrawingObj& RGBForFaceColor, 
                    const Color4DrawingObj& RGBForBoundaryColor,
                    const Color4DrawingObj& RGBForCenterPTColor);
    Geometry2DColor(const Geometry2DColor& geometry2DColor);

    Geometry2DColor& operator=(const Geometry2DColor& geometry2DColor);

    inline void set_RGB_for_face(const float& R, const float& G, const float& B, const float& alpha = 1.0);
    inline void set_RGB_for_boundary(const float& R, const float& G, const float& B, const float& alpha = 1.0);
    inline void set_RGB_for_center(const float& R, const float& G, const float& B, const float& alpha = 1.0);

    inline Color4DrawingObj get_RGB_for_face()     const { return m_RGBForFaceColor; };
    inline Color4DrawingObj get_RGB_for_boundary() const { return m_RGBForBoundaryColor; };
    inline Color4DrawingObj get_RGB_for_center()   const { return m_RGBForCenterPTColor; };

private:
    void copy_from(const Geometry2DColor& geometry2DColor);
};


inline void Geometry2DColor::set_RGB_for_face(const float& R, const float& G, const float& B, const float& alpha /*= 1.0*/)
{
    m_RGBForFaceColor.R = R;
    m_RGBForFaceColor.G = G;
    m_RGBForFaceColor.B = B;
    m_RGBForFaceColor.ALPHA = alpha;
}


inline void Geometry2DColor::set_RGB_for_boundary(const float& R, const float& G, const float& B, const float& alpha /*= 1.0*/)
{
    m_RGBForBoundaryColor.R = R;
    m_RGBForBoundaryColor.G = G;
    m_RGBForBoundaryColor.B = B;
    m_RGBForBoundaryColor.ALPHA = alpha;
}


inline void Geometry2DColor::set_RGB_for_center(const float& R, const float& G, const float& B, const float& alpha /*= 1.0*/)
{
    m_RGBForCenterPTColor.R = R;
    m_RGBForCenterPTColor.G = G;
    m_RGBForCenterPTColor.B = B;
    m_RGBForCenterPTColor.ALPHA = alpha;
}


#endif //_2D_GEOMETRY_COLOR_
