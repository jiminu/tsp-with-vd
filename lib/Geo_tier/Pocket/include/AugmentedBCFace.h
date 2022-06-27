#ifndef _AUGMENTEDBCFACE_H
#define _AUGMENTEDBCFACE_H

#include "BetaFace.h"

using namespace V::GeometryTier;

#include "AugmentedBCEdge.h"

//class AugmentedBCEdge;

class AugmentedBCFace : public BetaFace
{
private:
	rg_INT	m_boundingState;
	rg_BOOL m_isArtificial;

    rg_INT  m_indexOfExteriorRegionToLeftCell;
    rg_INT  m_indexOfExteriorRegionToRightCell;

public:
	AugmentedBCFace();
	AugmentedBCFace(const rg_INT& id);
	AugmentedBCFace(AugmentedBCEdge* e1, AugmentedBCEdge* e2, AugmentedBCEdge* e3 );
	AugmentedBCFace(const AugmentedBCFace& augBCFace);

	~AugmentedBCFace();


	rg_INT	getBoundingState() const;
	rg_BOOL isArtificial() const;
    rg_INT  getIndexOfExteriorRegionToLeftCell() const;
    rg_INT  getIndexOfExteriorRegionToRightCell() const;


	void setBoundingState(const rg_INT& boundingState);
	void setIsArtificial( const rg_BOOL& isArtificial);
    void setIndexOfExteriorRegionToLeftCell( const rg_INT& indexOfExteriorRegionToLeftCell );
    void setIndexOfExteriorRegionToRightCell( const rg_INT& indexOfExteriorRegionToRightCell );


	AugmentedBCFace& operator =(const AugmentedBCFace& augBCFace);
	
    void reverse();
	//void reverseEdgeOrientations();
	AugmentedBCEdge* findMateEdgeOfVertex(const AugmentedBCVertex* const vertex) const;
	AugmentedBCEdge* findIncidentEdgeOfVertexNotGivenEdge(const AugmentedBCVertex* const givenVertex, 
														  const AugmentedBCEdge* const notThisEdge) const;
};


#endif
