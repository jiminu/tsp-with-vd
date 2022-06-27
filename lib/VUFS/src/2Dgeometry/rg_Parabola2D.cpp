#include "rg_Parabola2D.h"



rg_Parabola2D::rg_Parabola2D()
{
}



rg_Parabola2D::rg_Parabola2D(const rg_Point2D& focus, const rg_Line2D& directrix)
: m_focus(focus)
, m_directrix(directrix)
{
}



rg_Parabola2D::rg_Parabola2D(const rg_Parabola2D& parabola)
: m_focus(parabola.m_focus)
, m_directrix(parabola.m_directrix)
{
}



rg_Parabola2D::~rg_Parabola2D()
{
}



rg_Parabola2D& rg_Parabola2D::operator =(const rg_Parabola2D& parabola)
{
    if (this != &parabola) {
        m_focus = parabola.m_focus;
        m_directrix = parabola.m_directrix;
    }

    return *this;
}



rg_Point2D  rg_Parabola2D::focus() const
{
    return m_focus;
}



rg_Line2D   rg_Parabola2D::directrix() const
{
    return m_directrix;
}



void        rg_Parabola2D::setFocus(const rg_Point2D& focus)
{
    m_focus = focus;
}



void        rg_Parabola2D::setDirectrix(const rg_Line2D& directrix)
{
    m_directrix = directrix;
}



void        rg_Parabola2D::setParabola(const rg_Point2D& focus, const rg_Line2D& directrix)
{
    m_focus = focus;
    m_directrix = directrix;
}



bool        rg_Parabola2D::convertIntoImplicitForm(double& A, double& B, double& C, double& D, double& E, double& F)
{
    bool isConversionComplete = true;

    if ( m_directrix.does_contain(m_focus) ) {
        isConversionComplete = false;
    }
    else {
        double a, b, c;
        m_directrix.convertIntoImplicitForm(a, b, c);

        double x_focus = m_focus.getX();
        double y_focus = m_focus.getY();

        A = a*a - 1;
        B = 2.0*a*b;
        C = b*b - 1;
        D = 2.0*a*c + 2.0*x_focus;
        E = 2.0*b*c + 2.0*y_focus;
        F = c*c - x_focus*x_focus - y_focus*y_focus;

        if ( rg_NEG(A) ) {
            A *= -1.0;
            B *= -1.0;
            C *= -1.0;
            D *= -1.0;
            E *= -1.0;
            F *= -1.0;
        }

        isConversionComplete = true;
    }

    return isConversionComplete;
}
