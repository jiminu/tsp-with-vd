#ifndef _AUGMENTEDBCEDGE_H
#define _AUGMENTEDBCEDGE_H

#include "BetaEdge.h"

using namespace V::GeometryTier;

class AugmentedBCVertex;
class AugmentedBCFace;
class AugmentedBCCell;

class AugmentedBCEdge : public BetaEdge
{
private:
	rg_INT	m_boundingState;
	rg_BOOL m_isArtificial;

    rg_INT  m_indexOfExteriorRegion;

public:
	AugmentedBCEdge();
	AugmentedBCEdge(const rg_INT& id);
	AugmentedBCEdge(AugmentedBCVertex* startVertex, AugmentedBCVertex* endVertex);
	AugmentedBCEdge(const AugmentedBCEdge& augBCVertex);

	~AugmentedBCEdge();


	rg_INT	getBoundingState() const;
	rg_BOOL isArtificial() const;
    rg_INT  getIndexOfExteriorRegion() const;


	void setBoundingState(const rg_INT& boundingState);
	void setIsArtificial( const rg_BOOL& isArtificial);
    void setIndexOfExteriorRegion( const rg_INT& indexOfExteriorRegion );

	AugmentedBCEdge& operator =(const AugmentedBCEdge& augBCVertex);

	AugmentedBCVertex* getBoundingVertexNotGivenVertex(const AugmentedBCVertex* const notThisVertex) const;
	
	AugmentedBCCell* getCCWNextCell(const AugmentedBCFace* const currFace) const;
    AugmentedBCFace* getCCWNextFace(const AugmentedBCFace* const currFace) const;
    AugmentedBCCell* getCWNextCell(const AugmentedBCFace* const currFace) const;
    AugmentedBCFace* getCWNextFace(const AugmentedBCFace* const currFace) const;
};

#endif
