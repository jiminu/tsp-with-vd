#ifndef _POINTERTOSIMPLEX_H
#define _POINTERTOSIMPLEX_H

#include "rg_Const.h"


const rg_INT VERTEX_SIMPLEX = 0;
const rg_INT EDGE_SIMPLEX   = 1;
const rg_INT FACE_SIMPLEX   = 2;
const rg_INT CELL_SIMPLEX	= 3;


class AugmentedBCVertex;
class AugmentedBCEdge;
class AugmentedBCFace;
class AugmentedBCCell;

namespace V {
    namespace GeometryTier {

class BetaVertex;
class BetaEdge;
class BetaFace;
class BetaCell;

    }
}

using namespace V::GeometryTier;

class PointerToSimplex
{
private:
    rg_INT  m_typeOfSimplex;
    void*   m_simplex;

public:
    PointerToSimplex();
    PointerToSimplex(const rg_INT& typeOfSimplex, void* simplex);
    PointerToSimplex(const PointerToSimplex& tempPointerToSimplex);
	PointerToSimplex(AugmentedBCVertex* vertex);
	PointerToSimplex(AugmentedBCEdge* edge);
	PointerToSimplex(AugmentedBCFace* face);
	PointerToSimplex(AugmentedBCCell* cell);
	PointerToSimplex(BetaVertex* vertex);
	PointerToSimplex(BetaEdge* edge);
	PointerToSimplex(BetaFace* face);
	PointerToSimplex(BetaCell* cell);
    ~PointerToSimplex();

    rg_INT      getTypeOfSimplex() const;
    void*       getSimplex() const;

    void        setTypeOfSimplex(const rg_INT& typeOfSimplex);
    void        setSimplex(void* simplex);
    void        setPointerToSimplex(const rg_INT& typeOfSimplex, void* simplex);

    PointerToSimplex& operator =(const PointerToSimplex& tempPointerToSimplex);
    PointerToSimplex& operator =(AugmentedBCVertex* vertex);
    PointerToSimplex& operator =(AugmentedBCEdge*   edge);
    PointerToSimplex& operator =(AugmentedBCFace*   face);
    PointerToSimplex& operator =(AugmentedBCCell*   cell);
};
#endif


