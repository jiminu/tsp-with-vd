#ifndef _AUGMENTEDBCVERTEX_H
#define _AUGMENTEDBCVERTEX_H

#include "BetaVertex.h"

using namespace V::GeometryTier;

class AugmentedBCVertex : public BetaVertex
{
private:
    BetaVertex* m_originalBetaVertex;

	rg_INT	    m_boundingState;
	rg_BOOL     m_isArtificial;

public:
	AugmentedBCVertex();
	AugmentedBCVertex(const rg_INT& id);
	AugmentedBCVertex(const AugmentedBCVertex& augBCVertex);

	~AugmentedBCVertex();

    BetaVertex* getOriginalBetaVertex() const;    
	rg_INT	    getBoundingState() const;
	rg_BOOL     isArtificial() const;

    void setOriginalBetaVertex(BetaVertex* originalBetaVertex);
	void setBoundingState(const rg_INT& boundingState);
	void setIsArtificial( const rg_BOOL& isArtificial);

	AugmentedBCVertex& operator =(const AugmentedBCVertex& augBCVertex);
    
};


#endif
