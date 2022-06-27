#ifndef _ROTAMER_H_
#define _ROTAMER_H_

#include "rg_Atom.h"
#include "Backbone.h"
#include "EnclosingSphereOfSpheres.h"
#include "Vector.h"

class V::GeometryTier::Residue;
typedef rg_Point3D Vector3D;

// This class represents a rotamer (rotational isomer).
// The actual dihedral angles of "m_residue" could be different from "m_dihedralAngles"
// because m_dihedralAngles are not reflected before those angles are fixed as final solution.

class Rotamer
{
private:
    V::GeometryTier::Residue* m_residue          ;
	rg_REAL* m_dihedralAngles   ;
	rg_INT   m_numDihedralAngles;
	inline void setInitialValues()
	{
		m_residue           = rg_NULL;	
		m_dihedralAngles    = rg_NULL;
		m_numDihedralAngles = 0;
	}
	void     setDihedralAngles(V::GeometryTier::Residue* residue);

public:
	Rotamer();
	Rotamer(const Rotamer& rotamer);
	Rotamer(V::GeometryTier::Residue* residue);
	~Rotamer();

	rg_REAL computeEnergyWithOtherRotamer( Rotamer& theOtherRotamer );
	rg_REAL computeEnergyWithBackbone    ( Backbone& backbone  );
	rg_REAL computeEnergyWithBackbone    ( Backbone* backbone  );

	rg_REAL computeEnergyWithOtherBackbone(V::GeometryTier::Residue* residue);
	rg_REAL computeEnergyWithOtherResidue     (V::GeometryTier::Residue* residue);
	rg_REAL computeEnergyWithItsBackbone      ();
	rg_REAL computeEnergyWithItsBackboneAtomsNOtherResidue     (V::GeometryTier::Residue* residue);

	// SCWRL3
	// The following functions are used for computing energy and based on the article
	// Yang Cao and Lin Song and Zhichao Miao and Yun Hu and Liqing Tian and and Taijiao Jiang: 
	// Improved side-chain modeling by coupling clash-detection guided iterative search with rotamer relaxation,
	// Bioinformatics, 27 (6), 785-790, 2011
	rg_REAL computeEnergyWithOtherRotamerBySCWRL3( Rotamer& theOtherRotamer );
	rg_REAL computeEnergyWithBackboneBySCWRL3    ( Backbone     & backbone     );
	rg_REAL computeRotamerPreferenceEnergyBySCWRL3(const rg_INDEX& rotLibID, const rg_INDEX& highestProbabilityRotLibID);

	void evaluateGradientVectorOfRotamerPreferenceEnergyForDihedralAnglesBySCWRL3( Vector& givenSolVector, const rg_INDEX& rotLibID, Vector& gradientVector );

	void computeMinimumEnclosingSphere(V::GeometryTier::EnclosingSphereOfSpheres& MES_R);
	void computeMinimumEnclosingSphereOfResidue(V::GeometryTier::EnclosingSphereOfSpheres& MES_R);
	rg_REAL computeSumOfPairwiseXVolumeWithOtherRotamer(Rotamer* rotamer);
	rg_REAL computeSumOfPairwiseXVolumeWithOtherRotamer(Rotamer* rotamer, rg_INT& numXAtoms);
	rg_REAL computeSumOfPairwiseXVolumeWithOtherBackbone(V::GeometryTier::Residue* residue, rg_INT& numXAtoms);
	rg_REAL computeSumOfPairwiseAtomXVolumesWithOtherResidue(V::GeometryTier::Residue* residue);
	rg_REAL computeSumOfPairwiseAtomXVolumesWithOtherResidue(V::GeometryTier::Residue* residue, rg_INT& numXAtoms);
	rg_REAL computeSumOfPairwiseAtomXVolumesWithItsBackbone();
	rg_REAL computeSumOfPairwiseAtomXVolumesWithItsBackbone(rg_INT& numXAtoms);
	rg_REAL computeSumOfPairwiseXVolumeWithItsBackboneAtomsNOtherResidue(V::GeometryTier::Residue* residue);
	rg_REAL computeSumOfPairwiseXVolumeWithProteinBackbone(Backbone& backbone);

	inline V::GeometryTier::Residue* getResidue() { return m_residue; }
	inline  rg_REAL* getDihedralAngles() { return m_dihedralAngles; }
	rg_INT  getDihedralAngles(rg_REAL*& dihedralAngles);
	inline  rg_INT  getNumDihedralAngles() const { return m_numDihedralAngles; }
	rg_INT  getAtomsOnSidechain(rg_dList<V::GeometryTier::Atom*>& atomsOnSidechain );
	rg_INT  getAtomsOnSidechain(V::GeometryTier::Atom**& atomsOnSidechain );
	rg_INT  getAtomsOfResidue( rg_dList<V::GeometryTier::Atom*>& atomsOfResidue );
	void    setDihedralAngles(const rg_REAL* dihedralAngles);
	void    setResidue(V::GeometryTier::Residue* residue);
	void    setResidueWithDihedralAngles(V::GeometryTier::Residue* residue);
	void    updateSidechainAtomCoordinates();
	void    updateSidechainAtomCoordinatesByApplyingRotLibID(const rg_INDEX& rotLibID);
	void    updateSidechainAtomCoordinatesByApplyingDihedralAngles(const rg_REAL* dihedralAngles);

    //void    perturbSidechainAtomsByRotation(const Vector3D& rotAxis, const rg_REAL& angle);

	Rotamer& operator=(const Rotamer& rotamer);

private:	
	void destroyDihedralAngles();
};

#endif