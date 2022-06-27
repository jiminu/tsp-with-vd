#include "AugmentedBCVertex.h"

AugmentedBCVertex::AugmentedBCVertex()
: BetaVertex()  
{
    m_originalBetaVertex = rg_NULL;

	m_boundingState = rg_UNKNOWN;	//EXTRANEOUS_SIMPLEX와 같으므로 주의.
	m_isArtificial  = rg_FALSE;
}


AugmentedBCVertex::AugmentedBCVertex(const rg_INT& id)
: BetaVertex(id)
{
    m_originalBetaVertex = rg_NULL;

	m_boundingState = rg_UNKNOWN;
	m_isArtificial  = rg_FALSE;
}


AugmentedBCVertex::AugmentedBCVertex(const AugmentedBCVertex& augBCVertex)
: BetaVertex( augBCVertex )
{
    m_originalBetaVertex = augBCVertex.m_originalBetaVertex;

	m_boundingState = augBCVertex.m_boundingState;
	m_isArtificial  = augBCVertex.m_isArtificial;
}


AugmentedBCVertex::~AugmentedBCVertex()
{
}


BetaVertex* AugmentedBCVertex::getOriginalBetaVertex() const
{
	return m_originalBetaVertex;
}


rg_INT	AugmentedBCVertex::getBoundingState() const
{
	return m_boundingState;
}

rg_BOOL AugmentedBCVertex::isArtificial() const
{
	return m_isArtificial;
}


void AugmentedBCVertex::setOriginalBetaVertex(BetaVertex* originalBetaVertex)
{
	m_originalBetaVertex = originalBetaVertex;
}

void AugmentedBCVertex::setBoundingState(const rg_INT& boundingState)
{
	m_boundingState = boundingState;
}


void AugmentedBCVertex::setIsArtificial( const rg_BOOL& isArtificial)
{
	m_isArtificial = isArtificial;
}


AugmentedBCVertex& AugmentedBCVertex::operator =(const AugmentedBCVertex& augBCVertex)
{
	if ( this == &augBCVertex )
		return *this;

	BetaVertex::operator =(augBCVertex);

    m_originalBetaVertex = augBCVertex.m_originalBetaVertex;

	m_boundingState = augBCVertex.m_boundingState;
	m_isArtificial  = augBCVertex.m_isArtificial;
	
	return *this;
}
