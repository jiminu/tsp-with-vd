#ifndef _NONMANIFOLDCORNER_H
#define _NONMANIFOLDCORNER_H

#include "rg_Const.h"
#include "PointerToSimplex.h"

const int EE_V_CONFIGURATION = 0;
const int EF_V_CONFIGURATION = 1;
const int ET_V_CONFIGURATION = 2;
const int FF_V_CONFIGURATION = 3;
const int FF_E_CONFIGURATION = 4;
const int FT_V_CONFIGURATION = 5;
const int FT_E_CONFIGURATION = 6;
const int TF_E_CONFIGURATION = 7;
const int TT_V_CONFIGURATION = 8;
const int TT_E_CONFIGURATION = 9;
const int TT_D_V_CONFIGURATION = 10;

class NonManifoldCorner
{
private:
    PointerToSimplex    m_NMSimplex;    //non-manifold simplex
    
    PointerToSimplex    m_firstSimplex;
    PointerToSimplex    m_secondSimplex;

    //rg_INT configurationStatus[3];

public:
    NonManifoldCorner();
    NonManifoldCorner(const PointerToSimplex& location,
                             const PointerToSimplex& firstSimplex, 
                             const PointerToSimplex& secondSimplex);
    NonManifoldCorner(AugmentedBCVertex* nm_simplex,
                             AugmentedBCEdge* firstSimplex, 
                             AugmentedBCEdge* secondSimplex);
    NonManifoldCorner(AugmentedBCVertex* nm_simplex,
                             AugmentedBCEdge* firstSimplex, 
                             AugmentedBCFace* secondSimplex);
    NonManifoldCorner(AugmentedBCVertex* nm_simplex,
                             AugmentedBCEdge* firstSimplex, 
                             AugmentedBCCell* secondSimplex);
    NonManifoldCorner(AugmentedBCVertex* nm_simplex,
                             AugmentedBCFace* firstSimplex, 
                             AugmentedBCFace* secondSimplex);
    NonManifoldCorner(AugmentedBCVertex* nm_simplex,
                             AugmentedBCFace* firstSimplex, 
                             AugmentedBCCell* secondSimplex);
    NonManifoldCorner(AugmentedBCVertex* nm_simplex,
                             AugmentedBCCell* firstSimplex, 
                             AugmentedBCCell* secondSimplex);
    NonManifoldCorner(AugmentedBCEdge* nm_simplex,
                             AugmentedBCFace* firstSimplex, 
                             AugmentedBCFace* secondSimplex);
    NonManifoldCorner(AugmentedBCEdge* location,
                             AugmentedBCFace* firstSimplex, 
                             AugmentedBCCell* secondSimplex);
    NonManifoldCorner(AugmentedBCEdge* nm_simplex,
                             AugmentedBCCell* firstSimplex, 
                             AugmentedBCFace* secondSimplex);
    NonManifoldCorner(AugmentedBCEdge* nm_simplex,
                             AugmentedBCCell* firstSimplex, 
                             AugmentedBCCell* secondSimplex);
    NonManifoldCorner(AugmentedBCVertex* nm_simplex);
    NonManifoldCorner(AugmentedBCFace* nm_simplex,
                             AugmentedBCCell* firstSimplex,
                             AugmentedBCCell* secondSimplex);
    
    NonManifoldCorner(const NonManifoldCorner& tempNonManifoldCorner);
    ~NonManifoldCorner();

    PointerToSimplex  getNMSimplex() const;
    PointerToSimplex* getpNMSimplex();      
    PointerToSimplex  getFirstSimplex() const;
    PointerToSimplex* getpFirstSimplex();
    PointerToSimplex  getSecondSimplex() const;
    PointerToSimplex* getpSecondSimplex();

    void*   getEdgeSimplex() const;
    void*   getFaceSimplex() const;
    void*   getCellSimplex() const;

    void setNMSimplex(const PointerToSimplex& location);
    void setFirstSimplex(const PointerToSimplex& firstSimplex);
    void setSecondSimplex(const PointerToSimplex& secondSimplex);
    void setSimplices(const PointerToSimplex& firstSimplex, const PointerToSimplex& secondSimplex);
    void setNonManifoldCorner(const PointerToSimplex& nm_simplex,
                                     const PointerToSimplex& firstSimplex, 
                                     const PointerToSimplex& secondSimplex);

    rg_INT getTypeOfNonManifoldConfiguration();

    NonManifoldCorner& operator =(const NonManifoldCorner& tempNonManifoldCorner);
};
#endif


