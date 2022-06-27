#ifndef _ANALYSISOFBETAUNIVERSE_H
#define _ANALYSISOFBETAUNIVERSE_H

#include "BetaUniverse.h"

#include <list>
using namespace std;


namespace V {

namespace GeometryTier {

class AnalysisOfBetaUniverse
{
private:
    BetaUniverse* m_betaUniverse;

public:
    AnalysisOfBetaUniverse();
    AnalysisOfBetaUniverse(BetaUniverse* betaUniverse);
    AnalysisOfBetaUniverse(const AnalysisOfBetaUniverse& analysis);
    ~AnalysisOfBetaUniverse();

    inline BetaUniverse* getBetaUniverse() { return m_betaUniverse; }
    inline void          setBetaUniverse(BetaUniverse* betaUniverse) { m_betaUniverse = betaUniverse; }

    AnalysisOfBetaUniverse& operator =(const AnalysisOfBetaUniverse& analysis);


    void reportQTStatistics(ofstream& fout);
    void reportBUStatistics(ofstream& fout, const rg_REAL& betaValue);
    void reportBUComparison(ofstream& fout, const rg_REAL& betaBefore, const rg_REAL& betaAfter);

    void reportGeometryInQT(ofstream& fout);

    //  for beta-span
    void    reportEdgesInBetaComplex(ofstream& fout, const rg_REAL& betaValue);


    //  for Quasi-triangulation
    void    countAvgNumOfVerticesEdgesFacesAndCellsIncidentToVertex(rg_REAL& avgNumNeighborVtxs, rg_REAL& avgNumIncidentEdges, rg_REAL& avgNumIncidentFaces, rg_REAL& avgNumIncidentCells) const;
    void    countAvgNumOfFacesAndCellsIncidentToEdge(rg_REAL& avgNumIncidentFaces, rg_REAL& avgNumIncidentCells) const;

    void    countCellsWithMinusVolumeInQT(rg_INT& numCellsWithMinusVolume) const;
    void    countMultiAdjacencyAnomaliesInQT(rg_INT& twoAdjacency, rg_INT& threeAdjacency, rg_INT& fourAdjacency) const;
    void    countGateEdgesInQT(rg_INT& numGateEdges) const;
    void    countDanglingFaceAnomaliesInQT(rg_INT& numDanglingFaceAnomalies) const;
    void    countInterWorldAnomaliesInQT(rg_INT& numInterWorldAnomalies) const;
    void    countCellsInWorldsInQT(rg_dList<rg_INT>& numCellsInWorlds) const;
    void    countDepthAndNumCellsOfWorldsInQT(rg_dList< pair<rg_INT, rg_INT> >& depthAndNumCellsOfWorlds) const;
    rg_INT  countMaximalDepthOfWorldsInQT() const;


    //  for beta-complex
    void    countVerticesInBetaComplex( const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const;
    void    countEdgesInBetaComplex(    const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const;
    void    countFacesInBetaComplex(    const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const;
    void    countCellsInBetaComplex(    const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numInterior) const;

    void    countCellsWithMinusVolumeInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numInterior) const;
    void    countCellsWithMinusVolumeAndNotMultiAdjacencyAnomalyInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numInterior, rg_REAL& sumInteriorCellVolume) const;
    void    countSingularEllipticFacesInBetaComplex(const rg_REAL& betaValue, rg_INT& numSingularEllipticFaces) const;

    void    countCellStateOf2AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numByTwoExtCells, rg_INT& numByExtNIntCells, rg_INT& numByTwoIntCells) const;
    void    countCellStateOf3AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numByTwoExtCells, rg_INT& numByExtNIntCells, rg_INT& numByTwoIntCells) const;
    void    countCellStateOf4AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numByTwoExtCells, rg_INT& numByExtNIntCells, rg_INT& numByTwoIntCells) const;

    void    countFacesIn2AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const;
    void    countFacesIn3AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const;
    void    countFacesIn4AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const;

    void    countGateEdgesInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const;
    void    countDanglingFaceAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const;

    void    countDanglingAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numDanglingFaceAnomaly, rg_INT& numDanglingCellAnomaly, rg_INT& numDanglingClusterAnomaly) const;

    void    countWorldsInBetaComplex(const rg_REAL& betaValue, rg_INT& numWorlds) const;


    void    countExtremeVtxBehaviorOnEdgesInBetaComplex( const rg_REAL& betaValue, rg_INT& numIrregularBehavior ) const;


    //  for comparison of 2 beta-complexes
    void    countVerticesWithChangedState( const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* numVerticesInChangedState) const; //  numVertices: int array with 16 elements.
    void    countEdgesWithChangedState(    const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* numEdgesInChangedState) const;    //  numEdges: int array with 16 elements.
    void    countFacesWithChangedState(    const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* numFacesInChangedState) const;    //  numFaces: int array with 16 elements.
    void    countCellsWithChangedState(    const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* numCellsInChangedState) const;    //  numCells: int array with 16 elements.

    void    count2AdjacencyAnomaliesWithChangedState( const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* num2AdjAnomalies) const; //  num2AdjAnomalies: int array with 9 elements.
    void    count3AdjacencyAnomaliesWithChangedState( const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* num3AdjAnomalies) const; //  num2AdjAnomalies: int array with 9 elements.
    void    count4AdjacencyAnomaliesWithChangedState( const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* num4AdjAnomalies) const; //  num2AdjAnomalies: int array with 9 elements.

    void    countGateEdgesWithChangedState(            const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* numGateEdges) const;
    void    countDanglingFaceAnomaliesWithChangedState(const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* numDanglingFaceAnomalies) const;
    void    countInterWorldAnomaliesWithChangedState(  const rg_REAL& beforeBetaValue, const rg_REAL& afterBetaValue, rg_INT* numInterWorldAnomalies) const;



    //  for geometry test in Quasi-triangulation
    void    evaluateAnglesOfCellsIncidentToEdgesInQT( list< pair<BetaEdge*, rg_REAL> >& angleSumOfCellsIncidentToEdge);
    void    checkDetailAroundEdgesInQT( ofstream& fout );


private:
    void    collectQFacesInWorlds(rg_dList<BetaFace*>& firstQFacesInWorlds) const;
};

} // namespace GeometryTier

} // namespace V

#endif
