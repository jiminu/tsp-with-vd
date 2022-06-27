#ifndef _AUGMENTEDBCCELL_H
#define _AUGMENTEDBCCELL_H

#include "BetaCell.h"

using namespace V::GeometryTier;

#include "AugmentedBCFace.h"
//class AugmentedBCFace;
class AugmentedBCEdge;
class AugmentedBCVertex;

class AugmentedBCCell : public BetaCell
{
private:
	rg_INT	m_boundingState;
	rg_BOOL m_isArtificial;

public:
	AugmentedBCCell();
	AugmentedBCCell(const rg_INT& id);
	AugmentedBCCell( AugmentedBCFace* f0, AugmentedBCFace* f1, AugmentedBCFace* f2, AugmentedBCFace* f3,   
					 AugmentedBCVertex* v0, AugmentedBCVertex* v1, AugmentedBCVertex* v2, AugmentedBCVertex* v3 );
	AugmentedBCCell( AugmentedBCFace* f0, AugmentedBCFace* f1, AugmentedBCFace* f2, AugmentedBCFace* f3 );
	AugmentedBCCell(const AugmentedBCCell& augBCCell);

	~AugmentedBCCell();


	rg_INT	getBoundingState() const;
	rg_BOOL isArtificial() const;


	void setBoundingState(const rg_INT& boundingState);
	void setIsArtificial( const rg_BOOL& isArtificial);

	AugmentedBCCell& operator =(const AugmentedBCCell& augBCCell);


	AugmentedBCFace* findNextIncidentFaceOfEdge(AugmentedBCEdge* givenEdge, AugmentedBCFace* firstIncidentFace);
	void findFacesIncidentToEdge(AugmentedBCEdge* givenEdge, AugmentedBCFace** incidentFaces);

    BetaVertex* findMateVertex( AugmentedBCFace* givenFace );
};


#endif
