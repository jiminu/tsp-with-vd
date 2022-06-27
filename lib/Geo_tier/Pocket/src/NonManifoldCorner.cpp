#include "NonManifoldCorner.h"

NonManifoldCorner::NonManifoldCorner()
{    
}



NonManifoldCorner::NonManifoldCorner(const PointerToSimplex& nmSimplex,
                         const PointerToSimplex& firstSimplex, const PointerToSimplex& secondSimplex)
{
    m_NMSimplex     = nmSimplex;
    m_firstSimplex  = firstSimplex;
    m_secondSimplex = secondSimplex;
}




NonManifoldCorner::NonManifoldCorner(const NonManifoldCorner& tempNonManifoldCorner)
{
    m_NMSimplex     = tempNonManifoldCorner.m_NMSimplex;
    m_firstSimplex  = tempNonManifoldCorner.m_firstSimplex;
    m_secondSimplex = tempNonManifoldCorner.m_secondSimplex;
}


NonManifoldCorner::NonManifoldCorner(AugmentedBCVertex* nmSimplex,
                                                     AugmentedBCEdge* firstSimplex, 
                                                     AugmentedBCEdge* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( VERTEX_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( EDGE_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( EDGE_SIMPLEX, (void*)secondSimplex );
}


NonManifoldCorner::NonManifoldCorner(AugmentedBCVertex* nmSimplex,
                                                     AugmentedBCEdge* firstSimplex, 
                                                     AugmentedBCFace* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( VERTEX_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( EDGE_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( FACE_SIMPLEX, (void*)secondSimplex );
}


NonManifoldCorner::NonManifoldCorner(AugmentedBCVertex* nmSimplex,
                                                     AugmentedBCEdge* firstSimplex, 
                                                     AugmentedBCCell* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( VERTEX_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( EDGE_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( CELL_SIMPLEX, (void*)secondSimplex );
}


NonManifoldCorner::NonManifoldCorner(AugmentedBCVertex* nmSimplex,
                                                     AugmentedBCFace* firstSimplex, 
                                                     AugmentedBCFace* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( VERTEX_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( FACE_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( FACE_SIMPLEX, (void*)secondSimplex );

}


NonManifoldCorner::NonManifoldCorner(AugmentedBCVertex* nmSimplex,
                                                     AugmentedBCFace* firstSimplex, 
                                                     AugmentedBCCell* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( VERTEX_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( FACE_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( CELL_SIMPLEX, (void*)secondSimplex );

}



NonManifoldCorner::NonManifoldCorner(AugmentedBCVertex* nmSimplex,
                                                     AugmentedBCCell* firstSimplex, 
                                                     AugmentedBCCell* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( VERTEX_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( CELL_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( CELL_SIMPLEX, (void*)secondSimplex );
}
NonManifoldCorner::NonManifoldCorner(AugmentedBCEdge* nmSimplex,
                                                     AugmentedBCFace* firstSimplex, 
                                                     AugmentedBCFace* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( EDGE_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( FACE_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( FACE_SIMPLEX, (void*)secondSimplex );
}



NonManifoldCorner::NonManifoldCorner(AugmentedBCEdge* nmSimplex,
                                                     AugmentedBCFace* firstSimplex, 
                                                     AugmentedBCCell* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( EDGE_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( FACE_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( CELL_SIMPLEX, (void*)secondSimplex );
}



NonManifoldCorner::NonManifoldCorner(AugmentedBCEdge* nmSimplex,
                                                     AugmentedBCCell* firstSimplex, 
                                                     AugmentedBCFace* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( EDGE_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( CELL_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( FACE_SIMPLEX, (void*)secondSimplex );
}



NonManifoldCorner::NonManifoldCorner(AugmentedBCEdge* nmSimplex,
                                                     AugmentedBCCell* firstSimplex, 
                                                     AugmentedBCCell* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( EDGE_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( CELL_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( CELL_SIMPLEX, (void*)secondSimplex );
}



NonManifoldCorner::NonManifoldCorner(AugmentedBCVertex* nmSimplex)
{
    m_NMSimplex     = PointerToSimplex( CELL_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( CELL_SIMPLEX, rg_NULL );
    m_secondSimplex = PointerToSimplex( CELL_SIMPLEX, rg_NULL );
}


NonManifoldCorner::NonManifoldCorner(AugmentedBCFace* nmSimplex,
                                                     AugmentedBCCell* firstSimplex, 
                                                     AugmentedBCCell* secondSimplex)
{
    m_NMSimplex     = PointerToSimplex( FACE_SIMPLEX, (void*)nmSimplex );
    m_firstSimplex  = PointerToSimplex( CELL_SIMPLEX, (void*)firstSimplex );
    m_secondSimplex = PointerToSimplex( CELL_SIMPLEX, (void*)secondSimplex );
}





NonManifoldCorner::~NonManifoldCorner()
{
}





PointerToSimplex  NonManifoldCorner::getNMSimplex() const
{
    return m_NMSimplex;
}

PointerToSimplex*  NonManifoldCorner::getpNMSimplex()
{
    return &m_NMSimplex;
}




PointerToSimplex NonManifoldCorner::getFirstSimplex() const
{
    return m_firstSimplex;
}

PointerToSimplex* NonManifoldCorner::getpFirstSimplex()
{
    return &m_firstSimplex;
}



PointerToSimplex NonManifoldCorner::getSecondSimplex() const
{
    return m_secondSimplex;
}

PointerToSimplex* NonManifoldCorner::getpSecondSimplex()
{
    return &m_secondSimplex;
}




void NonManifoldCorner::setNMSimplex(const PointerToSimplex& nmSimplex)
{
    m_NMSimplex = nmSimplex;
}




void NonManifoldCorner::setFirstSimplex(const PointerToSimplex& firstSimplex)
{
    m_firstSimplex  = firstSimplex;
}



void NonManifoldCorner::setSecondSimplex(const PointerToSimplex& secondSimplex)
{
    m_secondSimplex = secondSimplex;
}




void NonManifoldCorner::setSimplices(const PointerToSimplex& firstSimplex,
                                            const PointerToSimplex& secondSimplex)
{
    m_firstSimplex  = firstSimplex;
    m_secondSimplex = secondSimplex;
}




void NonManifoldCorner::setNonManifoldCorner(const PointerToSimplex& nmSimplex,
                         const PointerToSimplex& firstSimplex,
                         const PointerToSimplex& secondSimplex)
{
    m_NMSimplex     = nmSimplex;
    m_firstSimplex  = firstSimplex;
    m_secondSimplex = secondSimplex;
}



rg_INT NonManifoldCorner::getTypeOfNonManifoldConfiguration()
{
	if ( m_firstSimplex.getSimplex() == rg_NULL
		|| m_secondSimplex.getSimplex() == rg_NULL )
		return TT_D_V_CONFIGURATION;

    /*
    discriminationNumber is computed by 
        multiplying 3 numbers.
    3 numbers represent the simplex-classes of nm_simplex, 
                                               first simplex, 
                                               and second simplex, respectively.
    
    if simplex-class is vertex,      the number is 0+1(dimension of simplex + 1)
    if simplex-class is edge,        the number is 1+1
    if simplex-class is face,        the number is 2+1
    if simplex-class is tetrahedron, the number is 3+1

    discriminationNumber of configuration is summerized by
    EE_V configuration : (1+1) * (1+1) * (0+1) = 4
    EF_V configuration : (1+1) * (2+1) * (0+1) = 6
    ET_V configuration : (1+1) * (3+1) * (0+1) = 8
    FF_V configuration : (2+1) * (2+1) * (0+1) = 9
    FF_E configuration : (2+1) * (2+1) * (1+1) = 18
    FT_V configuration : (2+1) * (3+1) * (0+1) = 12
    FT_E configuration : (2+1) * (3+1) * (1+1) = 24
    TT_V configuration : (3+1) * (3+1) * (0+1) = 16
    TT_E configuration : (3+1) * (3+1) * (1+1) = 32
    TT_D_V configuration : (3+1) * (3+1) * (3+1) = 64 	
    */

    rg_INT discriminationNumber = ( m_NMSimplex.getTypeOfSimplex() + 1 )
                                * ( m_firstSimplex.getTypeOfSimplex() + 1 )
                                * ( m_secondSimplex.getTypeOfSimplex() + 1 );
    
    switch ( discriminationNumber )
    {
    case 4:     //EE_V configuration : (1+1) * (1+1) * (0+1) = 4
        return EE_V_CONFIGURATION;

    case 6:     //EF_V configuration : (1+1) * (2+1) * (0+1) = 6
        return EF_V_CONFIGURATION;

    case 8:     //ET_V configuration : (1+1) * (3+1) * (0+1) = 8
        return ET_V_CONFIGURATION;

    case 9:     //FF_V configuration : (2+1) * (2+1) * (0+1) = 9
        return FF_V_CONFIGURATION;

    case 12:    //FT_V configuration : (2+1) * (3+1) * (0+1) = 12
        return FT_V_CONFIGURATION;

    case 16:    //TT_V configuration : (3+1) * (3+1) * (0+1) = 16
        return TT_V_CONFIGURATION;

    case 18:    //FF_E configuration : (2+1) * (2+1) * (1+1) = 18
        return FF_E_CONFIGURATION;

    case 24:    //FT_E configuration : (2+1) * (3+1) * (1+1) = 24
        if ( m_firstSimplex.getTypeOfSimplex() == FACE_SIMPLEX )
            return FT_E_CONFIGURATION;
        else // if( m_firstSimplex.getTypeOfAugmentedBCSimplex() == CELL_SIMPLEX )
            return TF_E_CONFIGURATION;

    case 32:    //TT_E configuration : (3+1) * (3+1) * (1+1) = 32
        return TT_E_CONFIGURATION;

//     case 64:    //TT_D_V configuration : (3+1) * (3+1) * (3+1) = 64
//         return TT_D_V_CONFIGURATION;

    default:
        return -1;
    }
}




NonManifoldCorner& NonManifoldCorner::operator =(const NonManifoldCorner& tempNonManifoldCorner)
{
    if ( this == &tempNonManifoldCorner )
        return *this;

    m_NMSimplex     = tempNonManifoldCorner.m_NMSimplex;
    m_firstSimplex  = tempNonManifoldCorner.m_firstSimplex;
    m_secondSimplex = tempNonManifoldCorner.m_secondSimplex;

    return *this;
}






void* NonManifoldCorner::getEdgeSimplex() const
{
    if ( m_firstSimplex.getTypeOfSimplex() == EDGE_SIMPLEX ) {
        return m_firstSimplex.getSimplex();
    }
    else if ( m_secondSimplex.getTypeOfSimplex() == EDGE_SIMPLEX ) {
        return m_secondSimplex.getSimplex();
    }
    else {
        return rg_NULL;
    }

}



void* NonManifoldCorner::getFaceSimplex() const
{
    if ( m_firstSimplex.getTypeOfSimplex() == FACE_SIMPLEX ) {
        return m_firstSimplex.getSimplex();
    }
    else if ( m_secondSimplex.getTypeOfSimplex() == FACE_SIMPLEX ) {
        return m_secondSimplex.getSimplex();
    }
    else {
        return rg_NULL;
    }
}



void* NonManifoldCorner::getCellSimplex() const
{
    if ( m_firstSimplex.getTypeOfSimplex() == CELL_SIMPLEX ) {
        return m_firstSimplex.getSimplex();
    }
    else if ( m_secondSimplex.getTypeOfSimplex() == CELL_SIMPLEX ) {
        return m_secondSimplex.getSimplex();
    }
    else {
        return rg_NULL;
    }
}
