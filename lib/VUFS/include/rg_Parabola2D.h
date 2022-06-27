#ifndef RG_PARABOLA2D_H
#define RG_PARABOLA2D_H


#include "rg_Point2D.h"
#include "rg_Line2D.h"

class rg_Parabola2D
{
private:
    rg_Point2D  m_focus;
    rg_Line2D   m_directrix;

public:
    rg_Parabola2D();
    rg_Parabola2D(const rg_Point2D& focus, const rg_Line2D& directrix);
    rg_Parabola2D(const rg_Parabola2D& parabola);
    ~rg_Parabola2D();

    rg_Parabola2D& operator =(const rg_Parabola2D& parabola);

    rg_Point2D  focus() const;
    rg_Line2D   directrix() const;

    void        setFocus(const rg_Point2D& focus);
    void        setDirectrix(const rg_Line2D& directrix);
    void        setParabola(const rg_Point2D& focus, const rg_Line2D& directrix);

    /**
    * \brief convert a parabola  with a focus and a directrix into the implicit form ( \f$Ax^2 + Bxy + Cy^2 + Dx + Ey + C = 0 (A>=0)\f$ ).
    */
    bool        convertIntoImplicitForm(double& A, double& B, double& C, double& D, double& E, double& F);

};

#endif


