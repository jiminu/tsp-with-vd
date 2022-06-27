#include "AnalysisOfBetaUniverse.h"
using namespace V::GeometryTier;

#include <set>
using namespace std;


AnalysisOfBetaUniverse::AnalysisOfBetaUniverse()
: m_betaUniverse(rg_NULL)
{
}



AnalysisOfBetaUniverse::AnalysisOfBetaUniverse(BetaUniverse* betaUniverse)
: m_betaUniverse(betaUniverse)
{
}



AnalysisOfBetaUniverse::AnalysisOfBetaUniverse(const AnalysisOfBetaUniverse& analysis)
: m_betaUniverse(analysis.m_betaUniverse)
{
}



AnalysisOfBetaUniverse::~AnalysisOfBetaUniverse()
{
}





AnalysisOfBetaUniverse& AnalysisOfBetaUniverse::operator =(const AnalysisOfBetaUniverse& analysis)
{
    if ( this == &analysis ) {
        return *this;
    }

    m_betaUniverse = analysis.m_betaUniverse;

    return *this;
}





void AnalysisOfBetaUniverse::reportQTStatistics(ofstream& fout)
{
    rg_INT numCellsInQT = m_betaUniverse->getNumCells();
    rg_INT numFacesInQT = m_betaUniverse->getNumFaces();
    rg_INT numEdgesInQT = m_betaUniverse->getNumEdges();
    rg_INT numVerticesInQT = m_betaUniverse->getNumVertices();

    rg_REAL avgNumSimplexesIncidentToVertex[4] = {0.0, 0.0, 0.0, 0.0};
    rg_REAL avgNumSimplexesIncidentToEdge[2]   = {0.0, 0.0};

    //rg_INT numCellsWithMinusVolume = 0;
    rg_INT numMultiAdjacency[3] = { 0, 0, 0 };
    rg_INT numGateEdges = 0;
    rg_INT numDanglingFaceAnomaly = 0;
    rg_INT numInterWorldAnomaly = 0;


    countAvgNumOfVerticesEdgesFacesAndCellsIncidentToVertex( avgNumSimplexesIncidentToVertex[0], avgNumSimplexesIncidentToVertex[1], 
                                                             avgNumSimplexesIncidentToVertex[2], avgNumSimplexesIncidentToVertex[3]);
    countAvgNumOfFacesAndCellsIncidentToEdge(avgNumSimplexesIncidentToEdge[0], avgNumSimplexesIncidentToEdge[1]);

    //countCellsWithMinusVolumeInQT( numCellsWithMinusVolume );

    countMultiAdjacencyAnomaliesInQT( numMultiAdjacency[0], numMultiAdjacency[1], numMultiAdjacency[2] );
    countGateEdgesInQT( numGateEdges );
    countDanglingFaceAnomaliesInQT( numDanglingFaceAnomaly );
    countInterWorldAnomaliesInQT( numInterWorldAnomaly );

    rg_dList<rg_INT> numCellsInWorlds;
    countCellsInWorldsInQT( numCellsInWorlds );

    rg_INT maxDepthInQT = countMaximalDepthOfWorldsInQT();


    fout << numVerticesInQT << "\t";
    fout << numEdgesInQT << "\t";
    fout << numFacesInQT << "\t";
    fout << numCellsInQT << "\t";

    fout << avgNumSimplexesIncidentToVertex[0] << "\t";
    fout << avgNumSimplexesIncidentToVertex[1] << "\t";
    fout << avgNumSimplexesIncidentToVertex[2] << "\t";
    fout << avgNumSimplexesIncidentToVertex[3] << "\t";

    fout << avgNumSimplexesIncidentToEdge[0] << "\t";
    fout << avgNumSimplexesIncidentToEdge[1] << "\t";


    //fout << numCellsWithMinusVolume << "\t\t";    

    fout << numMultiAdjacency[0] << "\t";
    fout << numMultiAdjacency[1] << "\t";
    fout << numMultiAdjacency[2] << "\t";

    fout << numGateEdges << "\t";
    fout << numDanglingFaceAnomaly << "\t";
    fout << numCellsInWorlds.getSize() << "\t";
    fout << maxDepthInQT << "\t";

    numCellsInWorlds.reset4Loop();
    while ( numCellsInWorlds.setNext4Loop() ) {
        rg_INT numCellsOfThisWorld = numCellsInWorlds.getEntity();

        fout << numCellsOfThisWorld << "\t";
    }


    fout << endl;
}



void AnalysisOfBetaUniverse::reportBUStatistics(ofstream& fout, const rg_REAL& betaValue)
{
    rg_INT numSimplexes[4][4]        = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
    rg_INT numSingularEllipticFaces  = 0;
    rg_INT numMultiAdjacency[3][3]   = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
    rg_INT numFacesInMultiAdjacency[3][4]   = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
    rg_INT numGateEdges[4]           = {0, 0, 0, 0};
    rg_INT numDanglingFaceAnomaly[4] = {0, 0, 0, 0};
    rg_INT numInterWorldAnomaly[4]   = {0, 0, 0, 0};

    rg_INT numCellsWithMinusVolume[2] = {0, 0};
    rg_INT numCellsWithMinusVolumeAndNotAdjAnomaly[2] = {0, 0};
    rg_REAL sumInteriorCellVolume = 0.0;

    rg_INT numIrregularBehavior = 0;

    countVerticesInBetaComplex( betaValue, numSimplexes[0][0], numSimplexes[0][1], numSimplexes[0][2], numSimplexes[0][3]);
    countEdgesInBetaComplex(    betaValue, numSimplexes[1][0], numSimplexes[1][1], numSimplexes[1][2], numSimplexes[1][3]);
    countFacesInBetaComplex(    betaValue, numSimplexes[2][0], numSimplexes[2][1], numSimplexes[2][2], numSimplexes[2][3]);
    countCellsInBetaComplex(    betaValue, numSimplexes[3][0], numSimplexes[3][3]);

    countExtremeVtxBehaviorOnEdgesInBetaComplex( betaValue, numIrregularBehavior );

    countCellsWithMinusVolumeInBetaComplex( betaValue, numCellsWithMinusVolume[0], numCellsWithMinusVolume[1] );
    countCellsWithMinusVolumeAndNotMultiAdjacencyAnomalyInBetaComplex( betaValue, numCellsWithMinusVolumeAndNotAdjAnomaly[0], numCellsWithMinusVolumeAndNotAdjAnomaly[1], sumInteriorCellVolume );

    countSingularEllipticFacesInBetaComplex( betaValue, numSingularEllipticFaces );

    countCellStateOf2AdjacencyAnomaliesInBetaComplex( betaValue, numMultiAdjacency[0][0], numMultiAdjacency[0][1], numMultiAdjacency[0][2] );
    countCellStateOf3AdjacencyAnomaliesInBetaComplex( betaValue, numMultiAdjacency[1][0], numMultiAdjacency[1][1], numMultiAdjacency[1][2] );
    countCellStateOf4AdjacencyAnomaliesInBetaComplex( betaValue, numMultiAdjacency[2][0], numMultiAdjacency[2][1], numMultiAdjacency[2][2] );

    countFacesIn2AdjacencyAnomaliesInBetaComplex(betaValue, numFacesInMultiAdjacency[0][0], numFacesInMultiAdjacency[0][1], numFacesInMultiAdjacency[0][2], numFacesInMultiAdjacency[0][3]);
    countFacesIn3AdjacencyAnomaliesInBetaComplex(betaValue, numFacesInMultiAdjacency[1][0], numFacesInMultiAdjacency[1][1], numFacesInMultiAdjacency[1][2], numFacesInMultiAdjacency[1][3]);
    countFacesIn4AdjacencyAnomaliesInBetaComplex(betaValue, numFacesInMultiAdjacency[2][0], numFacesInMultiAdjacency[2][1], numFacesInMultiAdjacency[2][2], numFacesInMultiAdjacency[2][3]);

    countGateEdgesInBetaComplex(             betaValue, numGateEdges[0], numGateEdges[1], numGateEdges[2], numGateEdges[3] );
    countDanglingFaceAnomaliesInBetaComplex( betaValue, numDanglingFaceAnomaly[0], numDanglingFaceAnomaly[1], numDanglingFaceAnomaly[2], numDanglingFaceAnomaly[3] );

    rg_INT numWorlds;
    countWorldsInBetaComplex( betaValue, numWorlds );


    //fout << betaValue << "\t";

    // # exterior, singular, regular, and interior b-vertices
    fout << numSimplexes[0][0] << "\t" << numSimplexes[0][1] << "\t" << numSimplexes[0][2] << "\t" << numSimplexes[0][3] << "\t";

    // # exterior, singular, regular, and interior b-edges
    fout << numSimplexes[1][0] << "\t" << numSimplexes[1][1] << "\t" << numSimplexes[1][2] << "\t" << numSimplexes[1][3] << "\t";

    // # exterior, singular, regular, and interior b-faces
    fout << numSimplexes[2][0] << "\t" << numSimplexes[2][1] << "\t" << numSimplexes[2][2] << "\t" << numSimplexes[2][3] << "\t";

    // # exterior and interior b-cells
    fout << numSimplexes[3][0] << "\t" << numSimplexes[3][3] << "\t";


    // # 2-adjacency defined by two interior cell 
    fout << numMultiAdjacency[0][2] << "\t";
    
    // # 3-adjacency defined by two interior cell 
    fout << numMultiAdjacency[1][2] << "\t";
    
    // # 4-adjacency defined by two interior cell 
    fout << numMultiAdjacency[2][2] << "\t";
    

    // # gates : singular, regular, interior
    fout << numGateEdges[1] << "\t" << numGateEdges[2] << "\t" << numGateEdges[3] << "\t";

    // # dangling faces 
    fout << numDanglingFaceAnomaly[1] + numDanglingFaceAnomaly[2] + numDanglingFaceAnomaly[3] << "\t";

    // # world 
    fout << numWorlds << "\t";






    //fout << numIrregularBehavior << "\t\t";
    //fout << numSingularEllipticFaces << "\t\t";

    //fout << numCellsWithMinusVolume[0] << "\t" << numCellsWithMinusVolume[1] << "\t\t";
    //fout << numCellsWithMinusVolumeAndNotAdjAnomaly[0] << "\t" << numCellsWithMinusVolumeAndNotAdjAnomaly[1] << "\t" << sumInteriorCellVolume << "\t\t";

    //fout << numMultiAdjacency[0][0] << "\t" << numMultiAdjacency[0][1] << "\t" << numMultiAdjacency[0][2] << "\t\t";
    //fout << numFacesInMultiAdjacency[0][0] << "\t" << numFacesInMultiAdjacency[0][1] << "\t" << numFacesInMultiAdjacency[0][2] << "\t" << numFacesInMultiAdjacency[0][3] <<"\t\t";

    //fout << numMultiAdjacency[1][0] << "\t" << numMultiAdjacency[1][1] << "\t" << numMultiAdjacency[1][2] << "\t\t";
    //fout << numFacesInMultiAdjacency[1][0] << "\t" << numFacesInMultiAdjacency[1][1] << "\t" << numFacesInMultiAdjacency[1][2] << "\t" << numFacesInMultiAdjacency[1][3] <<"\t\t";

    //fout << numMultiAdjacency[2][0] << "\t" << numMultiAdjacency[2][1] << "\t" << numMultiAdjacency[2][2] << "\t\t";
    //fout << numFacesInMultiAdjacency[2][0] << "\t" << numFacesInMultiAdjacency[2][1] << "\t" << numFacesInMultiAdjacency[2][2] << "\t" << numFacesInMultiAdjacency[2][3] <<"\t\t";


    //fout << numGateEdges[0]           << "\t" << numGateEdges[1]           << "\t" << numGateEdges[2]           << "\t" << numGateEdges[3]           << "\t\t";
    //fout << numDanglingFaceAnomaly[0] << "\t" << numDanglingFaceAnomaly[1] << "\t" << numDanglingFaceAnomaly[2] << "\t" << numDanglingFaceAnomaly[3] << "\t\t";
    //fout << numInterWorldAnomaly[0]   << "\t" << numInterWorldAnomaly[1]   << "\t" << numInterWorldAnomaly[2]   << "\t" << numInterWorldAnomaly[3]   << "\t\t";

    fout << endl;
}



void AnalysisOfBetaUniverse::reportBUComparison(ofstream& fout, const rg_REAL& betaBefore, const rg_REAL& betaAfter)
{
    rg_INT numStateChangeInVertex[16] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    rg_INT numStateChangeInEdge[16]   = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    rg_INT numStateChangeInFace[16]   = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    rg_INT numStateChangeInCell[4]    = { 0, 0, 0, 0 };

    rg_INT numStateChangeIn2AdjacencyAnomaly[9]   = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };
    rg_INT numStateChangeIn3AdjacencyAnomaly[9]   = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };
    rg_INT numStateChangeIn4AdjacencyAnomaly[9]   = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

    countVerticesWithChangedState( betaBefore, betaAfter, numStateChangeInVertex); //  numVertices: int array with 16 elements.
    countEdgesWithChangedState(    betaBefore, betaAfter, numStateChangeInEdge);    //  numEdges: int array with 16 elements.
    countFacesWithChangedState(    betaBefore, betaAfter, numStateChangeInFace);    //  numFaces: int array with 16 elements.
    countCellsWithChangedState(    betaBefore, betaAfter, numStateChangeInCell);    //  numCells: int array with 16 elements.

    count2AdjacencyAnomaliesWithChangedState( betaBefore, betaAfter, numStateChangeIn2AdjacencyAnomaly); //  num2AdjAnomalies: int array with 9 elements.
    count3AdjacencyAnomaliesWithChangedState( betaBefore, betaAfter, numStateChangeIn3AdjacencyAnomaly); //  num2AdjAnomalies: int array with 9 elements.
    count4AdjacencyAnomaliesWithChangedState( betaBefore, betaAfter, numStateChangeIn4AdjacencyAnomaly); //  num2AdjAnomalies: int array with 9 elements.

    fout << "beta-value" << "\t" << betaBefore << "\t->\t" << betaAfter << endl;
    fout << "V" << "\t\t\t\t\t\t" 
         << "E" << "\t\t\t\t\t\t" 
         << "F" << "\t\t\t\t\t\t" 
         << "C" << "\t\t\t\t"
         << "2-adj." << "\t\t\t\t\t"
         << "3-adj." << "\t\t\t\t\t"
         << "4-adj." << endl;
    fout << "\t"<< "exterior" << "\t" << "singular" << "\t" << "regular" << "\t" << "interior" << "\t\t\t";
    fout << "exterior" << "\t" << "singular" << "\t" << "regular" << "\t" << "interior" << "\t\t\t";
    fout << "exterior" << "\t" << "singular" << "\t" << "regular" << "\t" << "interior" << "\t\t\t";
    fout << "exterior" << "\t" << "interior" << "\t\t\t";
    fout << "2-ext. cell" << "\t" << "ext-int cell" << "\t" << "2-int cell" << "\t\t\t";
    fout << "2-ext. cell" << "\t" << "ext-int cell" << "\t" << "2-int cell" << "\t\t\t";
    fout << "2-ext. cell" << "\t" << "ext-int cell" << "\t" << "2-int cell" << endl;

    rg_INT i=0;
    fout << "exterior" << "\t" ;
    for ( i=0; i<4; i++ ) {
        fout << numStateChangeInVertex[i] << "\t";
    }
    fout << "\t" << "exterior" << "\t" ;
    for ( i=0; i<4; i++ ) {
        fout << numStateChangeInEdge[i] << "\t";
    }
    fout << "\t" << "exterior" << "\t" ;
    for ( i=0; i<4; i++ ) {
        fout << numStateChangeInFace[i] << "\t";
    }
    fout << "\t" << "exterior" << "\t" ;
    for ( i=0; i<2; i++ ) {
        fout << numStateChangeInCell[i] << "\t";
    }

    fout << "\t" << "2-ext. cell" << "\t" ;
    for ( i=0; i<3; i++ ) {
        fout << numStateChangeIn2AdjacencyAnomaly[i] << "\t";
    }
    fout << "\t" << "2-ext. cell" << "\t" ;
    for ( i=0; i<3; i++ ) {
        fout << numStateChangeIn3AdjacencyAnomaly[i] << "\t";
    }
    fout << "\t" << "2-ext. cell" << "\t" ;
    for ( i=0; i<3; i++ ) {
        fout << numStateChangeIn4AdjacencyAnomaly[i] << "\t";
    }
    fout << endl;


    fout << "singular" << "\t" ;
    for ( i=4; i<8; i++ ) {
        fout << numStateChangeInVertex[i] << "\t";
    }
    fout << "\t" << "singular" << "\t" ;
    for ( i=4; i<8; i++ ) {
        fout << numStateChangeInEdge[i] << "\t";
    }
    fout << "\t" << "singular" << "\t" ;
    for ( i=4; i<8; i++ ) {
        fout << numStateChangeInFace[i] << "\t";
    }
    fout << "\t" << "interior" << "\t" ;
    for ( i=2; i<4; i++ ) {
        fout << numStateChangeInCell[i] << "\t";
    }

    fout << "\t" << "ext-int cell" << "\t" ;
    for ( i=3; i<6; i++ ) {
        fout << numStateChangeIn2AdjacencyAnomaly[i] << "\t";
    }
    fout << "\t" << "ext-int cell" << "\t" ;
    for ( i=3; i<6; i++ ) {
        fout << numStateChangeIn3AdjacencyAnomaly[i] << "\t";
    }
    fout << "\t" << "ext-int cell" << "\t" ;
    for ( i=3; i<6; i++ ) {
        fout << numStateChangeIn4AdjacencyAnomaly[i] << "\t";
    }
    fout << endl;


    fout << "regular" << "\t" ;
    for ( i=8; i<12; i++ ) {
        fout << numStateChangeInVertex[i] << "\t";
    }
    fout << "\t" << "regular" << "\t" ;
    for ( i=8; i<12; i++ ) {
        fout << numStateChangeInEdge[i] << "\t";
    }
    fout << "\t" << "regular" << "\t" ;
    for ( i=8; i<12; i++ ) {
        fout << numStateChangeInFace[i] << "\t";
    }
    fout << "\t\t";
    for ( i=0; i<2; i++ ) {
        fout << "\t";
    }

    fout << "\t" << "2-int cell" << "\t" ;
    for ( i=6; i<9; i++ ) {
        fout << numStateChangeIn2AdjacencyAnomaly[i] << "\t";
    }
    fout << "\t" << "2-int cell" << "\t" ;
    for ( i=6; i<9; i++ ) {
        fout << numStateChangeIn3AdjacencyAnomaly[i] << "\t";
    }
    fout << "\t" << "2-int cell" << "\t" ;
    for ( i=6; i<9; i++ ) {
        fout << numStateChangeIn4AdjacencyAnomaly[i] << "\t";
    }
    fout << endl;


    fout << "interior" << "\t" ;
    for ( i=12; i<16; i++ ) {
        fout << numStateChangeInVertex[i] << "\t";
    }
    fout << "\t" << "interior" << "\t" ;
    for ( i=12; i<16; i++ ) {
        fout << numStateChangeInEdge[i] << "\t";
    }
    fout << "\t" << "interior" << "\t" ;
    for ( i=12; i<16; i++ ) {
        fout << numStateChangeInFace[i] << "\t";
    }
    fout << endl << endl;
}



    
void AnalysisOfBetaUniverse::reportGeometryInQT(ofstream& fout)
{
    checkDetailAroundEdgesInQT( fout );


    //list< pair<BetaEdge*, rg_REAL> > angleSumOfCellsIncidentToEdge;
    //evaluateAnglesOfCellsIncidentToEdgesInQT( angleSumOfCellsIncidentToEdge );


    //list< pair<BetaEdge*, rg_REAL> >::const_iterator i_edge;
    //for ( i_edge = angleSumOfCellsIncidentToEdge.begin(); i_edge != angleSumOfCellsIncidentToEdge.end(); i_edge++ ) {
    //    BetaEdge* currEdge = i_edge->first;
    //    rg_REAL   angle    = i_edge->second;

    //    if ( currEdge->isOnConvexHull() ) {
    //        if ( rg_GT( angle, ANGLE_360_DEGREE ) ) {
    //            fout << currEdge->getID() << "\t" << "On SH" << "\t" <<  angle << endl;
    //        }
    //    }
    //    else {
    //        if ( rg_NE( angle, ANGLE_360_DEGREE ) ) {
    //            fout << currEdge->getID() << "\t" << "In SH" << "\t" <<  angle << endl;
    //        }
    //    }
    //}
}



//  for beta-span
void    AnalysisOfBetaUniverse::reportEdgesInBetaComplex(ofstream& fout, const rg_REAL& betaValue)
{
    rg_INT     numEdges  = m_betaUniverse->getNumEdges();
    BetaEdge** betaEdges = m_betaUniverse->getSortedEdges();

    for ( rg_INT i=0; i<numEdges; i++ ) {
        rg_INT boundingState = betaEdges[i]->getBoundingState(betaValue);

        //fout << "e" << betaEdges[i]->getID() << "\t";
        //fout << "v" << betaEdges[i]->getStartVertex()->getID() << "\t";
        //fout << "v" << betaEdges[i]->getEndVertex()->getID() << "\t";

        //fout << "a" << betaEdges[i]->getStartVertex()->getBallProperty()->getIDFromInput() << "\t";
        //fout << "a" << betaEdges[i]->getEndVertex()->getBallProperty()->getIDFromInput() << "\t";

        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                //fout << "exterior" << "\t";
                break;
            case SINGULAR_SIMPLEX:
                //fout << "singular" << "\t";
                break;
            case REGULAR_SIMPLEX:
                {
                fout << "e" << betaEdges[i]->getID() << "\t";
                fout << "v" << betaEdges[i]->getStartVertex()->getID() << "\t";
                fout << "v" << betaEdges[i]->getEndVertex()->getID() << "\t";

                fout << "a" << betaEdges[i]->getStartVertex()->getBallProperty()->getIDFromInput() << "\t";
                fout << "a" << betaEdges[i]->getEndVertex()->getBallProperty()->getIDFromInput() << "\t";

                fout << "regular" << "\t";

                const rg_BetaSpan* betaSpan = betaEdges[i]->getpBetaSpan();

                rg_INT   numInterval = betaSpan->getNumBetaInterval();
                rg_REAL* interval    = betaSpan->getBetaIntervals();
                for ( rg_INT i=0; i<=numInterval; i++ ) {
                    fout << interval[i] << "\t";
                }
                fout << endl;
                }
                break;
            case INTERIOR_SIMPLEX:
                //fout << "interior" << "\t";
                break;
        }

        //const BetaSpanCore* betaSpan = betaEdges[i]->getpBetaSpan();

        //rg_INT   numInterval = betaSpan->getNumBetaInterval();
        //rg_REAL* interval    = betaSpan->getBetaIntervals();
        //for ( rg_INT i=0; i<=numInterval; i++ ) {
        //    fout << interval[i] << "\t";
        //}
        //fout << endl;
    }
}



//  for Quasi-triangulation
void    AnalysisOfBetaUniverse::countCellsWithMinusVolumeInQT(rg_INT& numCellsWithMinusVolume) const
{
    numCellsWithMinusVolume = 0;

    rg_INT     numQTCell = m_betaUniverse->getNumCells();
    BetaCell** qtCell    = m_betaUniverse->getSortedCells();

    rg_INT i = 0;
    for ( i=0; i<numQTCell; i++ ) {
        rg_REAL signedVol = qtCell[i]->computeSignedVolume();

        if ( rg_LT( signedVol, 0.0 ) ) {
            numCellsWithMinusVolume++;
        }
    }
}



void    AnalysisOfBetaUniverse::countAvgNumOfVerticesEdgesFacesAndCellsIncidentToVertex(rg_REAL& avgNumNeighborVtxs, rg_REAL& avgNumIncidentEdges, rg_REAL& avgNumIncidentFaces, rg_REAL& avgNumIncidentCells) const
{
    rg_INT totalNumNeighborVtxs  = 0;
    rg_INT totalNumIncidentEdges = 0;
    rg_INT totalNumIncidentFaces = 0;
    rg_INT totalNumIncidentCells = 0;

    rg_INT       numQTVtx = m_betaUniverse->getNumVertices();
    BetaVertex** qtVertex = m_betaUniverse->getSortedVertices();
    for ( rg_INT i=0; i<numQTVtx; i++ ) {
        BetaVertex* currVtx = qtVertex[i];

        rg_dList<BetaVertex*> finiteVertexList;
        currVtx->searchFiniteVerticesInWholeWorld( finiteVertexList );

        rg_dList<BetaEdge*>   finiteEdgeList;
        currVtx->searchFiniteEdgesInWholeWorld(    finiteEdgeList   );

        rg_dList<BetaFace*>   finiteFaceList;
        currVtx->searchFiniteFacesInWholeWorld(    finiteFaceList   );

        rg_dList<BetaCell*>   finiteCellList;
        currVtx->searchFiniteCellsInWholeWorld(    finiteCellList   );

        totalNumNeighborVtxs  += finiteVertexList.getSize();
        totalNumIncidentEdges += finiteEdgeList.getSize();
        totalNumIncidentFaces += finiteFaceList.getSize();
        totalNumIncidentCells += finiteCellList.getSize();
    }

    avgNumNeighborVtxs  = totalNumNeighborVtxs/(double)numQTVtx;
    avgNumIncidentEdges = totalNumIncidentEdges/(double)numQTVtx;
    avgNumIncidentFaces = totalNumIncidentFaces/(double)numQTVtx;
    avgNumIncidentCells = totalNumIncidentCells/(double)numQTVtx;
}



void    AnalysisOfBetaUniverse::countAvgNumOfFacesAndCellsIncidentToEdge(rg_REAL& avgNumIncidentFaces, rg_REAL& avgNumIncidentCells) const
{
    rg_INT totalNumIncidentFaces = 0;
    rg_INT totalNumIncidentCells = 0;

    rg_INT       numQTEdges = m_betaUniverse->getNumEdges();
    BetaEdge**   qtEdge     = m_betaUniverse->getSortedEdges();
    for ( rg_INT i=0; i<numQTEdges; i++ ) {
        BetaEdge* currEdge = qtEdge[i];

        rg_dList<BetaFace*> finiteFaceList;
        currEdge->searchFiniteFacesInWholeWorld( finiteFaceList );

        rg_dList<BetaCell*> finiteCellList;
        currEdge->searchFiniteCellsInWholeWorld( finiteCellList );

        totalNumIncidentFaces += finiteFaceList.getSize();
        totalNumIncidentCells += finiteCellList.getSize();
    }

    avgNumIncidentFaces = totalNumIncidentFaces/(double)numQTEdges;
    avgNumIncidentCells = totalNumIncidentCells/(double)numQTEdges;
}



void    AnalysisOfBetaUniverse::countMultiAdjacencyAnomaliesInQT(rg_INT& twoAdjacency, rg_INT& threeAdjacency, rg_INT& fourAdjacency) const
{
    twoAdjacency   = 0;
    threeAdjacency = 0;
    fourAdjacency  = 0;

    rg_INT     numQTCell = m_betaUniverse->getNumCells();
    BetaCell** qtCell    = m_betaUniverse->getSortedCells();

    rg_INT i = 0;
    for ( i=0; i<numQTCell; i++ ) {
        if ( qtCell[i]->isVisited() ) {
            continue;
        }

        qtCell[i]->isVisited( rg_TRUE );

        rg_INT multiplicity = qtCell[i]->isMultiplicity();

        if ( multiplicity == 1 ) {
            continue;
        }
        else if ( multiplicity == 2 ) {
            twoAdjacency++;
            qtCell[i]->getNeighborCellWithMultiplicity()->isVisited( rg_TRUE );
        }
        else if ( multiplicity == 3 ) {
            threeAdjacency++;
            qtCell[i]->getNeighborCellWithMultiplicity()->isVisited( rg_TRUE );
        }
        else if ( multiplicity == 4 ) {
            fourAdjacency++;
            qtCell[i]->getNeighborCellWithMultiplicity()->isVisited( rg_TRUE );
        }
        else {
            // do nothing.
        }
    }

    for ( i=0; i<numQTCell; i++ ) {
        qtCell[i]->isVisited( rg_FALSE );
    }
}



void    AnalysisOfBetaUniverse::countGateEdgesInQT(rg_INT& numGateEdges) const
{
    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    rg_INT i = 0;
    for ( i=0; i<numQTEdge; i++ ) {
        if ( qtEdge[i]->isGateToSmallWorlds() ) {
            numGateEdges++;
        }
    }
}



void    AnalysisOfBetaUniverse::countDanglingFaceAnomaliesInQT(rg_INT& numDanglingFaceAnomalies) const
{
    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    rg_INT i = 0;
    for ( i=0; i<numQTEdge; i++ ) {
        if ( qtEdge[i]->isGateToSmallWorlds() ) {
            rg_dList<BetaFace*>* smallWorldList = qtEdge[i]->getSmallWorlds();

            BetaFace* currSmallWorld = rg_NULL;
            smallWorldList->reset4Loop();
            while ( smallWorldList->setNext4Loop() ) {
                currSmallWorld = smallWorldList->getEntity();

                if ( currSmallWorld->isIsolatedFace() ) {
                    numDanglingFaceAnomalies++;
                }
            }
        }
    }
}



void    AnalysisOfBetaUniverse::countInterWorldAnomaliesInQT(rg_INT& numInterWorldAnomalies) const
{
    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    rg_INT i = 0;
    for ( i=0; i<numQTEdge; i++ ) {
        if ( qtEdge[i]->isGateToSmallWorlds() ) {
            numInterWorldAnomalies += qtEdge[i]->getNumOfSmallWorlds();
        }
    }
}



void    AnalysisOfBetaUniverse::countCellsInWorldsInQT(rg_dList<rg_INT>& numCellsInWorlds) const
{
    rg_dList<BetaFace*> firstQFacesInWorlds;
    collectQFacesInWorlds( firstQFacesInWorlds );



    firstQFacesInWorlds.reset4Loop();
    while ( firstQFacesInWorlds.setNext4Loop() ) {
        BetaFace* seedFaceOfWorld = firstQFacesInWorlds.getEntity();

        if ( seedFaceOfWorld->isIsolatedFace() ) {
            numCellsInWorlds.add( 0 );
        }
        else {
            set<BetaCell*> cellsInThisWorld;

            rg_dList<BetaCell*> cellStack;
            cellStack.pushBack( seedFaceOfWorld->getLeftCell() );

	        while ( cellStack.getSize() > 0 ) {
		        BetaCell* popedCell = cellStack.popFront();
                cellsInThisWorld.insert( popedCell );

                BetaCell* neighbor[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
	            popedCell->getNeighborCells( neighbor );  
		        for ( rg_INT i=0; i<4; i++) {

                    if ( cellsInThisWorld.find( neighbor[i] ) == cellsInThisWorld.end() ) {
				        cellStack.pushBack( neighbor[i] );
                    }
		        }
            }


            rg_INT numCells = 0;
            set<BetaCell*>::iterator i_cell;
            for ( i_cell=cellsInThisWorld.begin(); i_cell!=cellsInThisWorld.end(); i_cell++ ) {
                BetaCell* currCell = *i_cell;
                if ( !currCell->isVirtual() ) {
                    numCells++;
                }
            }
            numCellsInWorlds.add( numCells );
        }
    }
}



void    AnalysisOfBetaUniverse::countDepthAndNumCellsOfWorldsInQT(rg_dList< pair<rg_INT, rg_INT> >& depthAndNumCellsOfWorlds) const
{
    rg_dList<BetaFace*> firstQFacesInWorlds;
    collectQFacesInWorlds( firstQFacesInWorlds );



    firstQFacesInWorlds.reset4Loop();
    while ( firstQFacesInWorlds.setNext4Loop() ) {
        BetaFace* seedFaceOfWorld = firstQFacesInWorlds.getEntity();

        rg_INT    depthOfThisWorld = seedFaceOfWorld->getDepth();

        if ( seedFaceOfWorld->isIsolatedFace() ) {
            depthAndNumCellsOfWorlds.add( make_pair(depthOfThisWorld, 0) );
        }
        else {
            set<BetaCell*> cellsInThisWorld;

            rg_dList<BetaCell*> cellStack;
            cellStack.pushBack( seedFaceOfWorld->getLeftCell() );

	        while ( cellStack.getSize() > 0 ) {
		        BetaCell* popedCell = cellStack.popFront();
                cellsInThisWorld.insert( popedCell );

                BetaCell* neighbor[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
	            popedCell->getNeighborCells( neighbor );  
		        for ( rg_INT i=0; i<4; i++) {

                    if ( cellsInThisWorld.find( neighbor[i] ) == cellsInThisWorld.end() ) {
				        cellStack.pushBack( neighbor[i] );
                    }
		        }
            }


            rg_INT numCells = 0;
            set<BetaCell*>::iterator i_cell;
            for ( i_cell=cellsInThisWorld.begin(); i_cell!=cellsInThisWorld.end(); i_cell++ ) {
                BetaCell* currCell = *i_cell;
                if ( !currCell->isVirtual() ) {
                    numCells++;
                }
            }
            depthAndNumCellsOfWorlds.add( make_pair(depthOfThisWorld, numCells) );
        }
    }
}


    
rg_INT    AnalysisOfBetaUniverse::countMaximalDepthOfWorldsInQT() const
{
    rg_INT maximalDepthOfWorld = 0;

    rg_dList<BetaFace*> firstQFacesInWorlds;
    collectQFacesInWorlds( firstQFacesInWorlds );
    
    firstQFacesInWorlds.reset4Loop();
    while ( firstQFacesInWorlds.setNext4Loop() ) {
        BetaFace* seedFaceOfWorld = firstQFacesInWorlds.getEntity();

        rg_INT    depthOfThisWorld = seedFaceOfWorld->getDepth();
        if ( maximalDepthOfWorld < depthOfThisWorld ) {
            maximalDepthOfWorld = depthOfThisWorld;
        }
    }

    return maximalDepthOfWorld;
}



//  for beta-complex
void    AnalysisOfBetaUniverse::countVerticesInBetaComplex( const rg_REAL& betaValue, 
                                                            rg_INT& numExterior, 
                                                            rg_INT& numSingular, 
                                                            rg_INT& numRegular, 
                                                            rg_INT& numInterior) const
{
    numExterior = 0;
    numSingular = 0;
    numRegular  = 0;
    numInterior = 0;

    rg_INT       numVertices  = m_betaUniverse->getNumVertices();
    BetaVertex** betaVertices = m_betaUniverse->getSortedVertices();

    for ( rg_INT i=0; i<numVertices; i++ ) {
        rg_INT boundingState = betaVertices[i]->getBoundingState(betaValue);

        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case SINGULAR_SIMPLEX:
                numSingular++;
                break;
            case REGULAR_SIMPLEX:
                numRegular++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                break;
        }
    }
}



void    AnalysisOfBetaUniverse::countEdgesInBetaComplex(    const rg_REAL& betaValue, 
                                                            rg_INT& numExterior, 
                                                            rg_INT& numSingular, 
                                                            rg_INT& numRegular, 
                                                            rg_INT& numInterior) const
{
    numExterior = 0;
    numSingular = 0;
    numRegular  = 0;
    numInterior = 0;

    rg_INT     numEdges  = m_betaUniverse->getNumEdges();
    BetaEdge** betaEdges = m_betaUniverse->getSortedEdges();

    for ( rg_INT i=0; i<numEdges; i++ ) {
        rg_INT boundingState = betaEdges[i]->getBoundingState(betaValue);

        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case SINGULAR_SIMPLEX:
                numSingular++;
                break;
            case REGULAR_SIMPLEX:
                numRegular++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                break;
        }
    }
}



void    AnalysisOfBetaUniverse::countFacesInBetaComplex(    const rg_REAL& betaValue, 
                                                            rg_INT& numExterior, 
                                                            rg_INT& numSingular, 
                                                            rg_INT& numRegular, 
                                                            rg_INT& numInterior) const
{
    numExterior = 0;
    numSingular = 0;
    numRegular  = 0;
    numInterior = 0;

    rg_INT     numFaces  = m_betaUniverse->getNumFaces();
    BetaFace** betaFaces = m_betaUniverse->getSortedFaces();

    for ( rg_INT i=0; i<numFaces; i++ ) {
        rg_INT boundingState = betaFaces[i]->getBoundingState(betaValue);

        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case SINGULAR_SIMPLEX:
                numSingular++;
                break;
            case REGULAR_SIMPLEX:
                numRegular++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                break;
        }
    }
}



void    AnalysisOfBetaUniverse::countCellsInBetaComplex(    const rg_REAL& betaValue, 
                                                            rg_INT& numExterior, 
                                                            rg_INT& numInterior) const
{
    numExterior = 0;
    numInterior = 0;

    rg_INT     numCells  = m_betaUniverse->getNumCells();
    BetaCell** betaCells = m_betaUniverse->getSortedCells();

    for ( rg_INT i=0; i<numCells; i++ ) {
        rg_INT boundingState = betaCells[i]->getBoundingState(betaValue);

        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                break;
        }
    }
}


    
void    AnalysisOfBetaUniverse::countCellsWithMinusVolumeInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numInterior) const
{
    numExterior = 0;
    numInterior = 0;

    rg_INT     numCells  = m_betaUniverse->getNumCells();
    BetaCell** betaCells = m_betaUniverse->getSortedCells();

    for ( rg_INT i=0; i<numCells; i++ ) {
        rg_REAL signedVol = betaCells[i]->computeSignedVolume();

        if ( rg_GE( signedVol, 0.0 ) ) {
            continue;
        }

        rg_INT boundingState = betaCells[i]->getBoundingState(betaValue);
        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                break;
        }
    }
}



void    AnalysisOfBetaUniverse::countCellsWithMinusVolumeAndNotMultiAdjacencyAnomalyInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numInterior, rg_REAL& sumInteriorCellVolume) const
{
    numExterior = 0;
    numInterior = 0;
    sumInteriorCellVolume = 0.0;

    rg_INT     numCells  = m_betaUniverse->getNumCells();
    BetaCell** betaCells = m_betaUniverse->getSortedCells();

    for ( rg_INT i=0; i<numCells; i++ ) {
        rg_REAL signedVol = betaCells[i]->computeSignedVolume();

        if ( rg_GE( signedVol, 0.0 ) ) {
            continue;
        }

        rg_INT multiplicity = betaCells[i]->isMultiplicity();
        if ( multiplicity != 1 ) {
            continue;
        }

        rg_INT boundingState = betaCells[i]->getBoundingState(betaValue);
        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                sumInteriorCellVolume += signedVol;
                break;
        }
    }
}



void    AnalysisOfBetaUniverse::countSingularEllipticFacesInBetaComplex(const rg_REAL& betaValue, rg_INT& numSingularEllipticFaces) const
{
    numSingularEllipticFaces = 0;

    rg_INT     numFaces  = m_betaUniverse->getNumFaces();
    BetaFace** betaFaces = m_betaUniverse->getSortedFaces();

    for ( rg_INT i=0; i<numFaces; i++ ) {
        if ( betaFaces[i]->isSingularEllipticFace(betaValue) ) {
            numSingularEllipticFaces++;
        }
    }
}



void    AnalysisOfBetaUniverse::countCellStateOf2AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, 
                                                                      rg_INT& numByTwoExtCells, 
                                                                      rg_INT& numByExtNIntCells, 
                                                                      rg_INT& numByTwoIntCells) const
{
    numByTwoExtCells  = 0;
    numByExtNIntCells = 0;
    numByTwoIntCells  = 0;

    rg_INT     numQTCell = m_betaUniverse->getNumCells();
    BetaCell** qtCell    = m_betaUniverse->getSortedCells();

    rg_INT i = 0;
    for ( i=0; i<numQTCell; i++ ) {
        if ( qtCell[i]->isVisited() ) {
            continue;
        }

        qtCell[i]->isVisited( rg_TRUE );

        rg_INT multiplicity = qtCell[i]->isMultiplicity();

        if ( multiplicity != 2 ) {
            continue;
        }

        BetaCell* neighbor = qtCell[i]->getNeighborCellWithMultiplicity();
        neighbor->isVisited( rg_TRUE );

        rg_INT stateOfCurr = qtCell[i]->getBoundingState( betaValue );
        rg_INT stateOfNeig = neighbor->getBoundingState( betaValue );

        if ( stateOfCurr == EXTERIOR_SIMPLEX && stateOfNeig == EXTERIOR_SIMPLEX ) {
            numByTwoExtCells++;
        }
        else if ( stateOfCurr == INTERIOR_SIMPLEX && stateOfNeig == INTERIOR_SIMPLEX ) {
            numByTwoIntCells++;
        }
        else {
            numByExtNIntCells++;
        }
    }

    for ( i=0; i<numQTCell; i++ ) {
        qtCell[i]->isVisited( rg_FALSE );
    }
}



void    AnalysisOfBetaUniverse::countCellStateOf3AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numByTwoExtCells, rg_INT& numByExtNIntCells, rg_INT& numByTwoIntCells) const
{
    numByTwoExtCells  = 0;
    numByExtNIntCells = 0;
    numByTwoIntCells  = 0;

    rg_INT     numQTCell = m_betaUniverse->getNumCells();
    BetaCell** qtCell    = m_betaUniverse->getSortedCells();

    rg_INT i = 0;
    for ( i=0; i<numQTCell; i++ ) {
        if ( qtCell[i]->isVisited() ) {
            continue;
        }

        qtCell[i]->isVisited( rg_TRUE );

        rg_INT multiplicity = qtCell[i]->isMultiplicity();

        if ( multiplicity != 3 ) {
            continue;
        }

        BetaCell* neighbor = qtCell[i]->getNeighborCellWithMultiplicity();
        neighbor->isVisited( rg_TRUE );

        rg_INT stateOfCurr = qtCell[i]->getBoundingState( betaValue );
        rg_INT stateOfNeig = neighbor->getBoundingState( betaValue );

        if ( stateOfCurr == EXTERIOR_SIMPLEX && stateOfNeig == EXTERIOR_SIMPLEX ) {
            numByTwoExtCells++;
        }
        else if ( stateOfCurr == INTERIOR_SIMPLEX && stateOfNeig == INTERIOR_SIMPLEX ) {
            numByTwoIntCells++;
        }
        else {
            numByExtNIntCells++;
        }
    }

    for ( i=0; i<numQTCell; i++ ) {
        qtCell[i]->isVisited( rg_FALSE );
    }
}



void    AnalysisOfBetaUniverse::countCellStateOf4AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numByTwoExtCells, rg_INT& numByExtNIntCells, rg_INT& numByTwoIntCells) const
{
    numByTwoExtCells  = 0;
    numByExtNIntCells = 0;
    numByTwoIntCells  = 0;

    rg_INT     numQTCell = m_betaUniverse->getNumCells();
    BetaCell** qtCell    = m_betaUniverse->getSortedCells();

    rg_INT i = 0;
    for ( i=0; i<numQTCell; i++ ) {
        if ( qtCell[i]->isVisited() ) {
            continue;
        }

        qtCell[i]->isVisited( rg_TRUE );

        rg_INT multiplicity = qtCell[i]->isMultiplicity();

        if ( multiplicity != 4 ) {
            continue;
        }

        BetaCell* neighbor = qtCell[i]->getNeighborCellWithMultiplicity();
        neighbor->isVisited( rg_TRUE );

        rg_INT stateOfCurr = qtCell[i]->getBoundingState( betaValue );
        rg_INT stateOfNeig = neighbor->getBoundingState( betaValue );

        if ( stateOfCurr == EXTERIOR_SIMPLEX && stateOfNeig == EXTERIOR_SIMPLEX ) {
            numByTwoExtCells++;
        }
        else if ( stateOfCurr == INTERIOR_SIMPLEX && stateOfNeig == INTERIOR_SIMPLEX ) {
            numByTwoIntCells++;
        }
        else {
            numByExtNIntCells++;
        }
    }

    for ( i=0; i<numQTCell; i++ ) {
        qtCell[i]->isVisited( rg_FALSE );
    }
}


void    AnalysisOfBetaUniverse::countFacesIn2AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, 
                                                                             rg_INT& numExterior, 
                                                                             rg_INT& numSingular, 
                                                                             rg_INT& numRegular, 
                                                                             rg_INT& numInterior) const
{
    numExterior = 0;
    numSingular = 0;
    numRegular  = 0;
    numInterior = 0;

    rg_INT     numFaces  = m_betaUniverse->getNumFaces();
    BetaFace** betaFaces = m_betaUniverse->getSortedFaces();

    for ( rg_INT i=0; i<numFaces; i++ ) {
        if ( !betaFaces[i]->isIn2AdjacencyAnomaly() ) {
            continue;
        }

        rg_INT boundingState = betaFaces[i]->getBoundingState(betaValue);

        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case SINGULAR_SIMPLEX:
                numSingular++;
                break;
            case REGULAR_SIMPLEX:
                numRegular++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                break;
        }
    }
}



void    AnalysisOfBetaUniverse::countFacesIn3AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, 
                                                                             rg_INT& numExterior, 
                                                                             rg_INT& numSingular, 
                                                                             rg_INT& numRegular, 
                                                                             rg_INT& numInterior) const
{
    numExterior = 0;
    numSingular = 0;
    numRegular  = 0;
    numInterior = 0;

    rg_INT     numFaces  = m_betaUniverse->getNumFaces();
    BetaFace** betaFaces = m_betaUniverse->getSortedFaces();

    for ( rg_INT i=0; i<numFaces; i++ ) {
        if ( !betaFaces[i]->isIn3AdjacencyAnomaly() ) {
            continue;
        }

        rg_INT boundingState = betaFaces[i]->getBoundingState(betaValue);

        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case SINGULAR_SIMPLEX:
                numSingular++;
                break;
            case REGULAR_SIMPLEX:
                numRegular++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                break;
        }
    }
}



void    AnalysisOfBetaUniverse::countFacesIn4AdjacencyAnomaliesInBetaComplex(const rg_REAL& betaValue, 
                                                                             rg_INT& numExterior, 
                                                                             rg_INT& numSingular, 
                                                                             rg_INT& numRegular, 
                                                                             rg_INT& numInterior) const
{
    numExterior = 0;
    numSingular = 0;
    numRegular  = 0;
    numInterior = 0;

    rg_INT     numFaces  = m_betaUniverse->getNumFaces();
    BetaFace** betaFaces = m_betaUniverse->getSortedFaces();

    for ( rg_INT i=0; i<numFaces; i++ ) {
        if ( !betaFaces[i]->isIn4AdjacencyAnomaly() ) {
            continue;
        }

        rg_INT boundingState = betaFaces[i]->getBoundingState(betaValue);

        switch ( boundingState ) {
            case EXTERIOR_SIMPLEX:
                numExterior++;
                break;
            case SINGULAR_SIMPLEX:
                numSingular++;
                break;
            case REGULAR_SIMPLEX:
                numRegular++;
                break;
            case INTERIOR_SIMPLEX:
                numInterior++;
                break;
        }
    }
}



void    AnalysisOfBetaUniverse::countGateEdgesInBetaComplex(const rg_REAL& betaValue, 
                                                            rg_INT& numExterior, 
                                                            rg_INT& numSingular, 
                                                            rg_INT& numRegular, 
                                                            rg_INT& numInterior) const
{
    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    rg_INT i = 0;
    for ( i=0; i<numQTEdge; i++ ) {
        if ( qtEdge[i]->isGateToSmallWorlds() ) {
            rg_INT boundingState = qtEdge[i]->getBoundingState(betaValue);

            switch ( boundingState ) {
                case EXTERIOR_SIMPLEX:
                    numExterior++;
                    break;
                case SINGULAR_SIMPLEX:
                    numSingular++;
                    break;
                case REGULAR_SIMPLEX:
                    numRegular++;
                    break;
                case INTERIOR_SIMPLEX:
                    numInterior++;
                    break;
            }        
        }
    }
}


    

void    AnalysisOfBetaUniverse::countDanglingFaceAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numExterior, rg_INT& numSingular, rg_INT& numRegular, rg_INT& numInterior) const
{
    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    rg_INT i = 0;
    for ( i=0; i<numQTEdge; i++ ) {
        if ( qtEdge[i]->isGateToSmallWorlds() ) {
            rg_dList<BetaFace*>* smallWorldList = qtEdge[i]->getSmallWorlds();

            BetaFace* currSmallWorld = rg_NULL;
            smallWorldList->reset4Loop();
            while ( smallWorldList->setNext4Loop() ) {
                currSmallWorld = smallWorldList->getEntity();

                if ( currSmallWorld->isIsolatedFace() ) {
                    rg_INT boundingState = currSmallWorld->getBoundingState(betaValue);

                    switch ( boundingState ) {
                        case EXTERIOR_SIMPLEX:
                            numExterior++;
                            break;
                        case SINGULAR_SIMPLEX:
                            numSingular++;
                            break;
                        case REGULAR_SIMPLEX:
                            numRegular++;
                            break;
                        case INTERIOR_SIMPLEX:
                            numInterior++;
                            break;
                    }        
                }
            }
        }
    }
}




void    AnalysisOfBetaUniverse::countDanglingAnomaliesInBetaComplex(const rg_REAL& betaValue, rg_INT& numDanglingFaceAnomaly, rg_INT& numDanglingCellAnomaly, rg_INT& numDanglingClusterAnomaly) const
{
    numDanglingFaceAnomaly = 0;
    numDanglingCellAnomaly = 0;
    numDanglingClusterAnomaly = 0;


    rg_dList<BetaFace*> firstQFacesInWorlds;
    collectQFacesInWorlds( firstQFacesInWorlds );



    firstQFacesInWorlds.reset4Loop();
    while ( firstQFacesInWorlds.setNext4Loop() ) {
        BetaFace* seedFaceOfWorld = firstQFacesInWorlds.getEntity();

        if ( seedFaceOfWorld->isVirtual() ) {
            // for root world.
            continue;
        }
        else if ( seedFaceOfWorld->isIsolatedFace() ) {
            if ( seedFaceOfWorld->getBoundingState( betaValue ) != EXTERIOR_SIMPLEX ) {
                numDanglingFaceAnomaly++;
            }
        }
        else {
            set<BetaCell*> cellsInThisWorld;

            rg_dList<BetaCell*> cellStack;
            cellStack.pushBack( seedFaceOfWorld->getLeftCell() );

	        while ( cellStack.getSize() > 0 ) {
		        BetaCell* popedCell = cellStack.popFront();
                cellsInThisWorld.insert( popedCell );

                BetaCell* neighbor[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
	            popedCell->getNeighborCells( neighbor );  
		        for ( rg_INT i=0; i<4; i++) {

                    if ( cellsInThisWorld.find( neighbor[i] ) == cellsInThisWorld.end() ) {
				        cellStack.pushBack( neighbor[i] );
                    }
		        }
            }


            rg_INT numInteriorCells = 0;
            set<BetaCell*>::iterator i_cell;
            for ( i_cell=cellsInThisWorld.begin(); i_cell!=cellsInThisWorld.end(); i_cell++ ) {
                BetaCell* currCell = *i_cell;

                if ( currCell->getBoundingState( betaValue ) == INTERIOR_SIMPLEX ) {
                    numInteriorCells++;
                }
            }

            if ( numInteriorCells > 0 ) {
                if ( cellsInThisWorld.size() == 2 && numInteriorCells == 2 ) {
                    numDanglingCellAnomaly++;
                }
                else if ( cellsInThisWorld.size() > 2 && numInteriorCells > 2 ) {
                    numDanglingClusterAnomaly++;
                }
            }
        }
    }
}



void    AnalysisOfBetaUniverse::countWorldsInBetaComplex(const rg_REAL& betaValue, rg_INT& numWorlds) const
{
    numWorlds = 0;


    rg_dList<BetaFace*> firstQFacesInWorlds;
    collectQFacesInWorlds( firstQFacesInWorlds );



    firstQFacesInWorlds.reset4Loop();
    while ( firstQFacesInWorlds.setNext4Loop() ) {
        BetaFace* seedFaceOfWorld = firstQFacesInWorlds.getEntity();

        if ( seedFaceOfWorld->isVirtual() ) {
            // for root world.
            numWorlds++;
        }
        else if ( seedFaceOfWorld->isIsolatedFace() ) {
            if ( seedFaceOfWorld->getBoundingState( betaValue ) != EXTERIOR_SIMPLEX ) {
                numWorlds++;
            }
        }
        else {
            set<BetaCell*> cellsInThisWorld;

            rg_dList<BetaCell*> cellStack;
            cellStack.pushBack( seedFaceOfWorld->getLeftCell() );

	        while ( cellStack.getSize() > 0 ) {
		        BetaCell* popedCell = cellStack.popFront();
                cellsInThisWorld.insert( popedCell );

                BetaCell* neighbor[4] = {rg_NULL, rg_NULL, rg_NULL, rg_NULL};
	            popedCell->getNeighborCells( neighbor );  
		        for ( rg_INT i=0; i<4; i++) {

                    if ( cellsInThisWorld.find( neighbor[i] ) == cellsInThisWorld.end() ) {
				        cellStack.pushBack( neighbor[i] );
                    }
		        }
            }

            set<BetaCell*>::iterator i_cell;
            for ( i_cell=cellsInThisWorld.begin(); i_cell!=cellsInThisWorld.end(); i_cell++ ) {
                BetaCell* currCell = *i_cell;

                if ( currCell->getBoundingState( betaValue ) != EXTERIOR_SIMPLEX ) {
                    numWorlds++;
                    break;
                }
            }
        }
    }
}



void    AnalysisOfBetaUniverse::countExtremeVtxBehaviorOnEdgesInBetaComplex( const rg_REAL& betaValue, rg_INT& numIrregularBehavior ) const
{
    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    rg_INT i = 0;
    for ( i=0; i<numQTEdge; i++ ) {
        BetaEdge* currEdge = qtEdge[i];

        rg_INT state = currEdge->getBoundingState( betaValue );
        if ( state == EXTERIOR_SIMPLEX || state == INTERIOR_SIMPLEX ) {
            continue;
        }

        Sphere  extremeBall[2] = { currEdge->getStartVertex()->getBall(), 
                                   currEdge->getEndVertex()->getBall() };
        extremeBall[0].resizeRadiusBy( betaValue );
        extremeBall[1].resizeRadiusBy( betaValue );

        rg_BOOL isIncludedInAnother[2] = { rg_FALSE, rg_FALSE };
        if ( extremeBall[1].doesContain( extremeBall[0].getCenter() ) ) {
            isIncludedInAnother[0] = rg_TRUE;
        }
        if ( extremeBall[0].doesContain( extremeBall[1].getCenter() ) ) {
            isIncludedInAnother[1] = rg_TRUE;
        }

        if ( isIncludedInAnother[0] || isIncludedInAnother[1] ) {
            numIrregularBehavior++;
        }
    }
}



//  for comparison of 2 beta-complex
//  numVertices: int array with 16 elements.
void    AnalysisOfBetaUniverse::countVerticesWithChangedState( const rg_REAL& beforeBetaValue, 
                                                               const rg_REAL& afterBetaValue, 
                                                               rg_INT* numVerticesInChangedState) const
{
    //  VInChange = numVerticesInChangedState
    //              after beta-value 
    //              exterior       singular       regular        interior
    //  exterior    VInChange[0]   VInChange[1]   VInChange[2]   VInChange[3]
    //  singular    VInChange[4]   VInChange[5]   VInChange[6]   VInChange[7]
    //  regular     VInChange[8]   VInChange[9]   VInChange[10]  VInChange[11]
    //  interior    VInChange[12]  VInChange[13]  VInChange[14]  VInChange[15]

    rg_INT i = 0;
    for ( i=0; i<16; i++ ) {
        numVerticesInChangedState[i] = 0;
    }

    rg_INT       numVertices  = m_betaUniverse->getNumVertices();
    BetaVertex** betaVertices = m_betaUniverse->getSortedVertices();

    for ( i=0; i<numVertices; i++ ) {
        rg_INT beforeState = betaVertices[i]->getBoundingState(beforeBetaValue);
        rg_INT afterState  = betaVertices[i]->getBoundingState(afterBetaValue);

        rg_INT i_before = 0;
        switch ( beforeState ) {
            case EXTERIOR_SIMPLEX:
                i_before = 0;
                break;
            case SINGULAR_SIMPLEX:
                i_before = 1;
                break;
            case REGULAR_SIMPLEX:
                i_before = 2;
                break;
            case INTERIOR_SIMPLEX:
                i_before = 3;
                break;
        }

        rg_INT i_after = 0;
        switch ( afterState ) {
            case EXTERIOR_SIMPLEX:
                i_after = 0;
                break;
            case SINGULAR_SIMPLEX:
                i_after = 1;
                break;
            case REGULAR_SIMPLEX:
                i_after = 2;
                break;
            case INTERIOR_SIMPLEX:
                i_after = 3;
                break;
        }

        rg_INT index = (4*i_before) + i_after;

        numVerticesInChangedState[ index ]++;
    }
}



//  numEdges: int array with 16 elements.
void    AnalysisOfBetaUniverse::countEdgesWithChangedState(    const rg_REAL& beforeBetaValue, 
                                                               const rg_REAL& afterBetaValue, 
                                                               rg_INT* numEdgesInChangedState) const
{
    //  EInChange = numEdgesInChangedState
    //              after beta-value 
    //              exterior       singular       regular        interior
    //  exterior    EInChange[0]   EInChange[1]   EInChange[2]   EInChange[3]
    //  singular    EInChange[4]   EInChange[5]   EInChange[6]   EInChange[7]
    //  regular     EInChange[8]   EInChange[9]   EInChange[10]  EInChange[11]
    //  interior    EInChange[12]  EInChange[13]  EInChange[14]  EInChange[15]

    rg_INT i = 0;
    for ( i=0; i<16; i++ ) {
        numEdgesInChangedState[i] = 0;
    }

    rg_INT       numEdges  = m_betaUniverse->getNumEdges();
    BetaEdge**   betaEdges = m_betaUniverse->getSortedEdges();

    for ( i=0; i<numEdges; i++ ) {
        rg_INT beforeState = betaEdges[i]->getBoundingState(beforeBetaValue);
        rg_INT afterState  = betaEdges[i]->getBoundingState(afterBetaValue);

        rg_INT i_before = 0;
        switch ( beforeState ) {
            case EXTERIOR_SIMPLEX:
                i_before = 0;
                break;
            case SINGULAR_SIMPLEX:
                i_before = 1;
                break;
            case REGULAR_SIMPLEX:
                i_before = 2;
                break;
            case INTERIOR_SIMPLEX:
                i_before = 3;
                break;
        }

        rg_INT i_after = 0;
        switch ( afterState ) {
            case EXTERIOR_SIMPLEX:
                i_after = 0;
                break;
            case SINGULAR_SIMPLEX:
                i_after = 1;
                break;
            case REGULAR_SIMPLEX:
                i_after = 2;
                break;
            case INTERIOR_SIMPLEX:
                i_after = 3;
                break;
        }


        if ( beforeState == REGULAR_SIMPLEX && afterState == INTERIOR_SIMPLEX ) {
            int stop = 1;
        }
        rg_INT index = (4*i_before) + i_after;

        numEdgesInChangedState[ index ]++;
    }
}



//  numFaces: int array with 16 elements.
void    AnalysisOfBetaUniverse::countFacesWithChangedState(    const rg_REAL& beforeBetaValue, 
                                                               const rg_REAL& afterBetaValue, 
                                                               rg_INT* numFacesInChangedState) const
{
    //  FInChange = numFacesInChangedState
    //              after beta-value 
    //              exterior       singular       regular        interior
    //  exterior    FInChange[0]   FInChange[1]   FInChange[2]   FInChange[3]
    //  singular    FInChange[4]   FInChange[5]   FInChange[6]   FInChange[7]
    //  regular     FInChange[8]   FInChange[9]   FInChange[10]  FInChange[11]
    //  interior    FInChange[12]  FInChange[13]  FInChange[14]  FInChange[15]

    rg_INT i = 0;
    for ( i=0; i<16; i++ ) {
        numFacesInChangedState[i] = 0;
    }

    rg_INT       numFaces  = m_betaUniverse->getNumFaces();
    BetaFace**   betaFaces = m_betaUniverse->getSortedFaces();

    for ( i=0; i<numFaces; i++ ) {
        rg_INT beforeState = betaFaces[i]->getBoundingState(beforeBetaValue);
        rg_INT afterState  = betaFaces[i]->getBoundingState(afterBetaValue);

        rg_INT i_before = 0;
        switch ( beforeState ) {
            case EXTERIOR_SIMPLEX:
                i_before = 0;
                break;
            case SINGULAR_SIMPLEX:
                i_before = 1;
                break;
            case REGULAR_SIMPLEX:
                i_before = 2;
                break;
            case INTERIOR_SIMPLEX:
                i_before = 3;
                break;
        }

        rg_INT i_after = 0;
        switch ( afterState ) {
            case EXTERIOR_SIMPLEX:
                i_after = 0;
                break;
            case SINGULAR_SIMPLEX:
                i_after = 1;
                break;
            case REGULAR_SIMPLEX:
                i_after = 2;
                break;
            case INTERIOR_SIMPLEX:
                i_after = 3;
                break;
        }

        rg_INT index = (4*i_before) + i_after;

        numFacesInChangedState[ index ]++;
    }
}



//  numCells: int array with 4 elements.
void    AnalysisOfBetaUniverse::countCellsWithChangedState(    const rg_REAL& beforeBetaValue, 
                                                               const rg_REAL& afterBetaValue, 
                                                               rg_INT* numCellsInChangedState) const    
{
    //  FInChange = numCellsInChangedState
    //              after beta-value 
    //              exterior       interior
    //  exterior    FInChange[0]   FInChange[1]
    //  interior    FInChange[2]   FInChange[3]

    rg_INT i = 0;
    for ( i=0; i<4; i++ ) {
        numCellsInChangedState[i] = 0;
    }

    rg_INT       numCells  = m_betaUniverse->getNumCells();
    BetaCell**   betaCells = m_betaUniverse->getSortedCells();

    for ( i=0; i<numCells; i++ ) {
        rg_INT beforeState = betaCells[i]->getBoundingState(beforeBetaValue);
        rg_INT afterState  = betaCells[i]->getBoundingState(afterBetaValue);

        rg_INT i_before = 0;
        switch ( beforeState ) {
            case EXTERIOR_SIMPLEX:
                i_before = 0;
                break;
            case INTERIOR_SIMPLEX:
                i_before = 1;
                break;
        }

        rg_INT i_after = 0;
        switch ( afterState ) {
            case EXTERIOR_SIMPLEX:
                i_after = 0;
                break;
            case INTERIOR_SIMPLEX:
                i_after = 1;
                break;
        }

        rg_INT index = (2*i_before) + i_after;

        numCellsInChangedState[ index ]++;
    }
}



//  num2AdjAnomalies: int array with 9 elements.
void    AnalysisOfBetaUniverse::count2AdjacencyAnomaliesWithChangedState( const rg_REAL& beforeBetaValue, 
                                                                          const rg_REAL& afterBetaValue, 
                                                                          rg_INT* num2AdjAnomalies) const
{
    //                      after beta-value 
    //                      2 ext. cells          ext. & int. cells     2 int. cells        
    //  2 ext. cells        num2AdjAnomalies[0]   num2AdjAnomalies[1]   num2AdjAnomalies[2]   
    //  ext. & int. cells   num2AdjAnomalies[3]   num2AdjAnomalies[4]   num2AdjAnomalies[5]   
    //  2 int. cells        num2AdjAnomalies[6]   num2AdjAnomalies[7]   num2AdjAnomalies[8]  

    rg_INT i = 0;
    for ( i=0; i<9; i++ ) {
        num2AdjAnomalies[i] = 0;
    }

    rg_INT     numQTCell = m_betaUniverse->getNumCells();
    BetaCell** qtCell    = m_betaUniverse->getSortedCells();

    for ( i=0; i<numQTCell; i++ ) {
        BetaCell* currCell = qtCell[i];
        if ( currCell->isVisited() ) {
            continue;
        }

        currCell->isVisited( rg_TRUE );

        rg_INT multiplicity = currCell->isMultiplicity();

        if ( multiplicity != 2 ) {
            continue;
        }

        BetaCell* neighbor = currCell->getNeighborCellWithMultiplicity();
        neighbor->isVisited( rg_TRUE );

        rg_INT beforeStateOfCurr = currCell->getBoundingState( beforeBetaValue );
        rg_INT beforeStateOfNeig = neighbor->getBoundingState( beforeBetaValue );

        rg_INT i_before = 0;

        if ( beforeStateOfCurr == EXTERIOR_SIMPLEX && beforeStateOfNeig == EXTERIOR_SIMPLEX ) {
            i_before = 0;
        }
        else if ( beforeStateOfCurr == INTERIOR_SIMPLEX && beforeStateOfNeig == INTERIOR_SIMPLEX ) {
            i_before = 2;
        }
        else {
            i_before = 1;
        }


        rg_INT afterStateOfCurr = currCell->getBoundingState( afterBetaValue );
        rg_INT afterStateOfNeig = neighbor->getBoundingState( afterBetaValue );

        rg_INT i_after = 0;

        if ( afterStateOfCurr == EXTERIOR_SIMPLEX && afterStateOfNeig == EXTERIOR_SIMPLEX ) {
            i_after = 0;
        }
        else if ( afterStateOfCurr == INTERIOR_SIMPLEX && afterStateOfNeig == INTERIOR_SIMPLEX ) {
            i_after = 2;
        }
        else {
            i_after = 1;
        }

        rg_INT index = (3*i_before) + i_after;
        num2AdjAnomalies[ index ]++;
    }

    for ( i=0; i<numQTCell; i++ ) {
        qtCell[i]->isVisited( rg_FALSE );
    }
}



//  num2AdjAnomalies: int array with 9 elements.
void    AnalysisOfBetaUniverse::count3AdjacencyAnomaliesWithChangedState( const rg_REAL& beforeBetaValue, 
                                                                          const rg_REAL& afterBetaValue, 
                                                                          rg_INT* num3AdjAnomalies) const
{
    //                      after beta-value 
    //                      2 ext. cells          ext. & int. cells     2 int. cells        
    //  2 ext. cells        num3AdjAnomalies[0]   num3AdjAnomalies[1]   num3AdjAnomalies[2]   
    //  ext. & int. cells   num3AdjAnomalies[3]   num3AdjAnomalies[4]   num3AdjAnomalies[5]   
    //  2 int. cells        num3AdjAnomalies[6]   num3AdjAnomalies[7]   num3AdjAnomalies[8]  

    rg_INT i = 0;
    for ( i=0; i<9; i++ ) {
        num3AdjAnomalies[i] = 0;
    }

    rg_INT     numQTCell = m_betaUniverse->getNumCells();
    BetaCell** qtCell    = m_betaUniverse->getSortedCells();

    for ( i=0; i<numQTCell; i++ ) {
        BetaCell* currCell = qtCell[i];
        if ( currCell->isVisited() ) {
            continue;
        }

        currCell->isVisited( rg_TRUE );

        rg_INT multiplicity = currCell->isMultiplicity();

        if ( multiplicity != 3 ) {
            continue;
        }

        BetaCell* neighbor = currCell->getNeighborCellWithMultiplicity();
        neighbor->isVisited( rg_TRUE );

        rg_INT beforeStateOfCurr = currCell->getBoundingState( beforeBetaValue );
        rg_INT beforeStateOfNeig = neighbor->getBoundingState( beforeBetaValue );

        rg_INT i_before = 0;

        if ( beforeStateOfCurr == EXTERIOR_SIMPLEX && beforeStateOfNeig == EXTERIOR_SIMPLEX ) {
            i_before = 0;
        }
        else if ( beforeStateOfCurr == INTERIOR_SIMPLEX && beforeStateOfNeig == INTERIOR_SIMPLEX ) {
            i_before = 2;
        }
        else {
            i_before = 1;
        }


        rg_INT afterStateOfCurr = currCell->getBoundingState( afterBetaValue );
        rg_INT afterStateOfNeig = neighbor->getBoundingState( afterBetaValue );

        rg_INT i_after = 0;

        if ( afterStateOfCurr == EXTERIOR_SIMPLEX && afterStateOfNeig == EXTERIOR_SIMPLEX ) {
            i_after = 0;
        }
        else if ( afterStateOfCurr == INTERIOR_SIMPLEX && afterStateOfNeig == INTERIOR_SIMPLEX ) {
            i_after = 2;
        }
        else {
            i_after = 1;
        }

        rg_INT index = (3*i_before) + i_after;
        num3AdjAnomalies[ index ]++;
    }

    for ( i=0; i<numQTCell; i++ ) {
        qtCell[i]->isVisited( rg_FALSE );
    }
}



//  num2AdjAnomalies: int array with 9 elements.
void    AnalysisOfBetaUniverse::count4AdjacencyAnomaliesWithChangedState( const rg_REAL& beforeBetaValue, 
                                                                          const rg_REAL& afterBetaValue, 
                                                                          rg_INT* num4AdjAnomalies) const
{
    //                      after beta-value 
    //                      2 ext. cells          ext. & int. cells     2 int. cells        
    //  2 ext. cells        num4AdjAnomalies[0]   num4AdjAnomalies[1]   num4AdjAnomalies[2]   
    //  ext. & int. cells   num4AdjAnomalies[3]   num4AdjAnomalies[4]   num4AdjAnomalies[5]   
    //  2 int. cells        num4AdjAnomalies[6]   num4AdjAnomalies[7]   num4AdjAnomalies[8]  

    rg_INT i = 0;
    for ( i=0; i<9; i++ ) {
        num4AdjAnomalies[i] = 0;
    }

    rg_INT     numQTCell = m_betaUniverse->getNumCells();
    BetaCell** qtCell    = m_betaUniverse->getSortedCells();

    for ( i=0; i<numQTCell; i++ ) {
        BetaCell* currCell = qtCell[i];
        if ( currCell->isVisited() ) {
            continue;
        }

        currCell->isVisited( rg_TRUE );

        rg_INT multiplicity = currCell->isMultiplicity();

        if ( multiplicity != 4 ) {
            continue;
        }

        BetaCell* neighbor = currCell->getNeighborCellWithMultiplicity();
        neighbor->isVisited( rg_TRUE );

        rg_INT beforeStateOfCurr = currCell->getBoundingState( beforeBetaValue );
        rg_INT beforeStateOfNeig = neighbor->getBoundingState( beforeBetaValue );

        rg_INT i_before = 0;

        if ( beforeStateOfCurr == EXTERIOR_SIMPLEX && beforeStateOfNeig == EXTERIOR_SIMPLEX ) {
            i_before = 0;
        }
        else if ( beforeStateOfCurr == INTERIOR_SIMPLEX && beforeStateOfNeig == INTERIOR_SIMPLEX ) {
            i_before = 2;
        }
        else {
            i_before = 1;
        }


        rg_INT afterStateOfCurr = currCell->getBoundingState( afterBetaValue );
        rg_INT afterStateOfNeig = neighbor->getBoundingState( afterBetaValue );

        rg_INT i_after = 0;

        if ( afterStateOfCurr == EXTERIOR_SIMPLEX && afterStateOfNeig == EXTERIOR_SIMPLEX ) {
            i_after = 0;
        }
        else if ( afterStateOfCurr == INTERIOR_SIMPLEX && afterStateOfNeig == INTERIOR_SIMPLEX ) {
            i_after = 2;
        }
        else {
            i_after = 1;
        }

        rg_INT index = (3*i_before) + i_after;
        num4AdjAnomalies[ index ]++;
    }

    for ( i=0; i<numQTCell; i++ ) {
        qtCell[i]->isVisited( rg_FALSE );
    }
}



void    AnalysisOfBetaUniverse::countGateEdgesWithChangedState(            const rg_REAL& beforeBetaValue, 
                                                                           const rg_REAL& afterBetaValue, 
                                                                           rg_INT* numGateEdges) const
{
    //              after beta-value 
    //              exterior       singular       regular        interior
    //  exterior    numGateEdges[0]   numGateEdges[1]   numGateEdges[2]   numGateEdges[3]
    //  singular    numGateEdges[4]   numGateEdges[5]   numGateEdges[6]   numGateEdges[7]
    //  regular     numGateEdges[8]   numGateEdges[9]   numGateEdges[10]  numGateEdges[11]
    //  interior    numGateEdges[12]  numGateEdges[13]  numGateEdges[14]  numGateEdges[15]

    rg_INT i = 0;
    for ( i=0; i<16; i++ ) {
        numGateEdges[i] = 0;
    }

    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    for ( i=0; i<numQTEdge; i++ ) {
        if ( qtEdge[i]->isGateToSmallWorlds() ) {
            rg_INT beforeState = qtEdge[i]->getBoundingState(beforeBetaValue);
            rg_INT afterState  = qtEdge[i]->getBoundingState(afterBetaValue);

            rg_INT i_before = 0;
            switch ( beforeState ) {
                case EXTERIOR_SIMPLEX:
                    i_before = 0;
                    break;
                case SINGULAR_SIMPLEX:
                    i_before = 1;
                    break;
                case REGULAR_SIMPLEX:
                    i_before = 2;
                    break;
                case INTERIOR_SIMPLEX:
                    i_before = 3;
                    break;
            }                    
            
            rg_INT i_after = 0;
            switch ( afterState ) {
                case EXTERIOR_SIMPLEX:
                    i_after = 0;
                    break;
                case SINGULAR_SIMPLEX:
                    i_after = 1;
                    break;
                case REGULAR_SIMPLEX:
                    i_after = 2;
                    break;
                case INTERIOR_SIMPLEX:
                    i_after = 3;
                    break;
            }        

            rg_INT index = (4*i_before) + i_after;

            numGateEdges[ index ]++;
        }
    }
}



void    AnalysisOfBetaUniverse::countDanglingFaceAnomaliesWithChangedState(const rg_REAL& beforeBetaValue, 
                                                                           const rg_REAL& afterBetaValue, 
                                                                           rg_INT* numDanglingFaceAnomalies) const
{
    //  N_DFA = numDanglingFaceAnomalies
    //              after beta-value 
    //              exterior   singular   regular    interior
    //  exterior    N_DFA[0]   N_DFA[1]   N_DFA[2]   N_DFA[3]
    //  singular    N_DFA[4]   N_DFA[5]   N_DFA[6]   N_DFA[7]
    //  regular     N_DFA[8]   N_DFA[9]   N_DFA[10]  N_DFA[11]
    //  interior    N_DFA[12]  N_DFA[13]  N_DFA[14]  N_DFA[15]

    rg_INT i = 0;
    for ( i=0; i<16; i++ ) {
        numDanglingFaceAnomalies[i] = 0;
    }


    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    for ( i=0; i<numQTEdge; i++ ) {
        if ( qtEdge[i]->isGateToSmallWorlds() ) {
            rg_dList<BetaFace*>* smallWorldList = qtEdge[i]->getSmallWorlds();

            BetaFace* currSmallWorld = rg_NULL;
            smallWorldList->reset4Loop();
            while ( smallWorldList->setNext4Loop() ) {
                currSmallWorld = smallWorldList->getEntity();

                if ( currSmallWorld->isIsolatedFace() ) {
                    rg_INT beforeState = currSmallWorld->getBoundingState(beforeBetaValue);
                    rg_INT afterState  = currSmallWorld->getBoundingState(afterBetaValue);

                    rg_INT i_before = 0;
                    switch ( beforeState ) {
                        case EXTERIOR_SIMPLEX:
                            i_before = 0;
                            break;
                        case SINGULAR_SIMPLEX:
                            i_before = 1;
                            break;
                        case REGULAR_SIMPLEX:
                            i_before = 2;
                            break;
                        case INTERIOR_SIMPLEX:
                            i_before = 3;
                            break;
                    }                    
                    
                    rg_INT i_after = 0;
                    switch ( afterState ) {
                        case EXTERIOR_SIMPLEX:
                            i_after = 0;
                            break;
                        case SINGULAR_SIMPLEX:
                            i_after = 1;
                            break;
                        case REGULAR_SIMPLEX:
                            i_after = 2;
                            break;
                        case INTERIOR_SIMPLEX:
                            i_after = 3;
                            break;
                    }        

                    rg_INT index = (4*i_before) + i_after;

                    numDanglingFaceAnomalies[ index ]++;
                }
            }
        }
    }
}



void    AnalysisOfBetaUniverse::countInterWorldAnomaliesWithChangedState(  const rg_REAL& beforeBetaValue, 
                                                                           const rg_REAL& afterBetaValue, 
                                                                           rg_INT* numInterWorldAnomalies) const
{
}



//  for geometry test in Quasi-triangulation
void    AnalysisOfBetaUniverse::evaluateAnglesOfCellsIncidentToEdgesInQT( list< pair<BetaEdge*, rg_REAL> >& angleSumOfCellsIncidentToEdge)
{
    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    for ( rg_INT i=0; i<numQTEdge; i++ ) {
        BetaEdge* currEdge = qtEdge[i];

        rg_Point3D  startPoint = currEdge->getStartVertex()->getBall().getCenter();
        rg_Point3D  endPoint   = currEdge->getEndVertex()->getBall().getCenter();
        rg_Point3D  normal     = endPoint - startPoint;
        normal.normalize();

        Plane       basePlane(normal, startPoint);


        rg_REAL angleByIncidentCells = 0.0;

        rg_dList<BetaCell*> finiteCellList;
        currEdge->searchFiniteCellsInWholeWorld( finiteCellList );

        finiteCellList.reset4Loop();
        while ( finiteCellList.setNext4Loop() ) {
            BetaCell* currCell = finiteCellList.getEntity();

            BetaFace* incidentFace[2] = {rg_NULL, rg_NULL};
            currCell->findFacesIncidentToEdge(currEdge, incidentFace);

	        BetaVertex* mateVertex[2] = {incidentFace[0]->findMateVertexOfEdge(currEdge), incidentFace[1]->findMateVertexOfEdge(currEdge)};
	        rg_Point3D  matePoint[2]  = {mateVertex[0]->getBall().getCenter(),        mateVertex[1]->getBall().getCenter()};
            rg_Point3D  projectedMatePoint[2]  = { basePlane.projectPointOnPlane(matePoint[0]), basePlane.projectPointOnPlane(matePoint[1])};

            rg_REAL angle = basePlane.computeAngleInCCW(projectedMatePoint[0], startPoint, projectedMatePoint[1]);

            if ( currCell->computeSignedVolume() >= 0.0 ) {
	            angleByIncidentCells += angle;
            }
            else {
                angleByIncidentCells -= (ANGLE_360_DEGREE-angle);
            }
        }


        angleSumOfCellsIncidentToEdge.push_back( make_pair( currEdge, angleByIncidentCells ) );
    }
}



void    AnalysisOfBetaUniverse::checkDetailAroundEdgesInQT( ofstream& fout )
{
    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    for ( rg_INT i=0; i<numQTEdge; i++ ) {
        BetaEdge* currEdge = qtEdge[i];


        rg_Point3D  startPoint = currEdge->getStartVertex()->getBall().getCenter();
        rg_Point3D  endPoint   = currEdge->getEndVertex()->getBall().getCenter();
        rg_Point3D  normal     = endPoint - startPoint;
        normal.normalize();
        Plane       basePlane(normal, startPoint);


        rg_REAL angleByIncidentCells = 0.0;

        rg_dList<BetaCell*> finiteCellList;
        currEdge->searchFiniteCellsInWholeWorld( finiteCellList );


        finiteCellList.reset4Loop();
        while ( finiteCellList.setNext4Loop() ) {
            BetaCell* currCell = finiteCellList.getEntity();

            BetaFace* incidentFace[2] = {rg_NULL, rg_NULL};
            currCell->findFacesIncidentToEdge(currEdge, incidentFace);

	        BetaVertex* mateVertex[2] = {incidentFace[0]->findMateVertexOfEdge(currEdge), incidentFace[1]->findMateVertexOfEdge(currEdge)};
	        rg_Point3D  matePoint[2]  = {mateVertex[0]->getBall().getCenter(),        mateVertex[1]->getBall().getCenter()};
            rg_Point3D  projectedMatePoint[2]  = { basePlane.projectPointOnPlane(matePoint[0]), basePlane.projectPointOnPlane(matePoint[1])};

            rg_REAL angle = basePlane.computeAngleInCCW(projectedMatePoint[0], startPoint, projectedMatePoint[1]);

            rg_REAL volume = currCell->computeSignedVolume();
            if ( rg_LT( volume, 0.0 ) ) {
                rg_REAL betaValueIntoInterior = currCell->getMinTangentSphere().getRadius();

                rg_BOOL isCommonFace[2] = {rg_FALSE, rg_FALSE};
                for ( rg_INT j=0; j<2; j++ ) {
                    Interval_REAL regularInterval;
                    incidentFace[j]->getpBetaSpan()->findBetaIntervalOfGivenState(REGULAR_SIMPLEX, regularInterval);

                    if ( betaValueIntoInterior == regularInterval.getLowerValue() ) {
                        isCommonFace[j] = rg_FALSE;
                    }
                    else if ( betaValueIntoInterior == regularInterval.getUpperValue() ) {
                        isCommonFace[j] = rg_TRUE;
                    }
                }

                if ( isCommonFace[0] != isCommonFace[1] ) {
                    angle = -(ANGLE_360_DEGREE - angle);
                }
            }

            angleByIncidentCells += angle;

            //angleByIncidentCells += angle;

            //if ( currCell->computeSignedVolume() >= 0.0 ) {
	           // angleByIncidentCells += angle;
            //}
            //else {
            //    angleByIncidentCells -= (ANGLE_360_DEGREE-angle);
            //}
        }

        if ( currEdge->isOnConvexHull() ) {
                fout << "e_" << currEdge->getID() << "\t" << "On SH" << "\t" <<  angleByIncidentCells << endl;
            //if ( rg_GT( angleByIncidentCells, ANGLE_360_DEGREE ) ) {
            //    fout << "e_" << currEdge->getID() << "\t" << "On SH" << "\t" <<  angleByIncidentCells << endl;
            //}
            //else {
            //    continue;
            //}
        }
        else {
                fout << "e_" << currEdge->getID() << "\t" << "In SH" << "\t" <<  angleByIncidentCells << endl;
            //if ( rg_NE( angleByIncidentCells, ANGLE_360_DEGREE ) ) {
            //    fout << "e_" << currEdge->getID() << "\t" << "In SH" << "\t" <<  angleByIncidentCells << endl;
            //}
            //else {
            //    continue;
            //}
        }

        // continue;
        Ball* startPtrBall = currEdge->getStartVertex()->getBallProperty();
        Ball* endPtrBall   = currEdge->getEndVertex()->getBallProperty();

        Sphere startBall = currEdge->getStartVertex()->getBall();
        Sphere endBall   = currEdge->getEndVertex()->getBall();
        fout << "\t" << "v_" << currEdge->getStartVertex()->getID() 
             << "(b_" << startPtrBall->getIDFromInput() << ")"
             << "\t" << startBall.getCenter().getX() << "\t" << startBall.getCenter().getY() << "\t" << startBall.getCenter().getZ() 
             << "\t" << startBall.getRadius() << endl;
        fout << "\t" << "v_" << currEdge->getEndVertex()->getID() 
             << "(b_" << endPtrBall->getIDFromInput() << ")"
             << "\t" << endBall.getCenter().getX() << "\t" << endBall.getCenter().getY() << "\t" << endBall.getCenter().getZ() 
             << "\t" << endBall.getRadius() << endl;
        fout << endl;



        rg_dList<BetaFace*> finiteFaceList;
        currEdge->searchFiniteFacesInWholeWorld( finiteFaceList );
        
        fout << "\t" << finiteFaceList.getSize() << " faces" << endl;        
        finiteFaceList.reset4Loop();
        while ( finiteFaceList.setNext4Loop() ) {
            BetaFace*   currFace = finiteFaceList.getEntity();

            BetaVertex* mateVtx  = currFace->findMateVertexOfEdge( currEdge );
            Sphere      mateBall = mateVtx->getBall();


            fout << "\t" << "f_" << currFace->getID();
            if ( currFace->getEdgeOrientation( currEdge ) ) {
                fout << "\t" << "(+)";
            }
            else {
                fout << "\t" << "(-)";
            }

            if ( currFace->getLeftCell()->isVirtual() ) {
                fout << "\t" << "c_" << currFace->getLeftCell()->getID() << "*";
            }
            else {
                fout << "\t" << "c_" << currFace->getLeftCell()->getID();
            }
            if ( currFace->getRightCell()->isVirtual() ) {
                fout << "\t" << "c_" << currFace->getRightCell()->getID() << "*";
            }
            else {
                fout << "\t" << "c_" << currFace->getRightCell()->getID();
            }

            fout << "\t" << "v_" << mateVtx->getID() 
                 << "(b_" << mateVtx->getBallProperty()->getIDFromInput() << ")"
                 << "\t" << mateBall.getCenter().getX() << "\t" << mateBall.getCenter().getY() << "\t" << mateBall.getCenter().getZ() 
                 << "\t" << mateBall.getRadius();

            rg_Point3D  prjectedMatePoint = basePlane.projectPointOnPlane( mateBall.getCenter() );
            fout << "\t" << prjectedMatePoint.getX() << "\t" << prjectedMatePoint.getY() << "\t" << prjectedMatePoint.getZ();
            fout << endl;


        }
        fout << endl;


        fout << "\t" << finiteCellList.getSize() << " cells" << endl;
        finiteCellList.reset4Loop();
        while ( finiteCellList.setNext4Loop() ) {
            BetaCell* currCell = finiteCellList.getEntity();

            fout << "\t" << "c_" << currCell->getID();

            BetaFace* incidentFace[2] = {rg_NULL, rg_NULL};
            currCell->findFacesIncidentToEdge(currEdge, incidentFace);

	        BetaVertex* mateVertex[2] = {incidentFace[0]->findMateVertexOfEdge(currEdge), incidentFace[1]->findMateVertexOfEdge(currEdge)};
	        rg_Point3D  matePoint[2]  = {mateVertex[0]->getBall().getCenter(),        mateVertex[1]->getBall().getCenter()};
            rg_Point3D  projectedMatePoint[2]  = { basePlane.projectPointOnPlane(matePoint[0]), basePlane.projectPointOnPlane(matePoint[1])};

            rg_REAL angle = basePlane.computeAngleInCCW(projectedMatePoint[0], startPoint, projectedMatePoint[1]);

            rg_REAL volume = currCell->computeSignedVolume();
            if ( rg_LT( volume, 0.0 ) ) {
                rg_REAL betaValueIntoInterior = currCell->getMinTangentSphere().getRadius();

                rg_BOOL isCommonFace[2] = {rg_FALSE, rg_FALSE};
                for ( rg_INT j=0; j<2; j++ ) {
                    Interval_REAL regularInterval;
                    incidentFace[j]->getpBetaSpan()->findBetaIntervalOfGivenState(REGULAR_SIMPLEX, regularInterval);

                    if ( betaValueIntoInterior == regularInterval.getLowerValue() ) {
                        isCommonFace[j] = rg_FALSE;
                    }
                    else if ( betaValueIntoInterior == regularInterval.getUpperValue() ) {
                        isCommonFace[j] = rg_TRUE;
                    }
                }

                if ( isCommonFace[0] != isCommonFace[1] ) {
                    angle = -(ANGLE_360_DEGREE - angle);
                }
            }

            //if ( currCell->computeSignedVolume() < 0.0 ) {
            //    angle = -(ANGLE_360_DEGREE-angle);
            //}

            fout << "\t" << "f_" << incidentFace[0]->getID();
            fout << "\t" << "f_" << incidentFace[1]->getID();
            fout << "\t" << angle;

        	BetaVertex** vtx = currCell->getVertices();

            fout << "\t" << "v_" << vtx[0]->getID();
            fout << "\t" << "v_" << vtx[1]->getID();
            fout << "\t" << "v_" << vtx[2]->getID();
            fout << "\t" << "v_" << vtx[3]->getID();

            fout << "\t" << currCell->computeSignedVolume();

            Sphere tgnSphere = currCell->getMinTangentSphere();
            fout << "\t" << tgnSphere.getCenter().getX() << "\t" << tgnSphere.getCenter().getY() << "\t" << tgnSphere.getCenter().getZ() 
                 << "\t" << tgnSphere.getRadius();

            fout << endl;
        }


        fout << endl << endl;
    }
}
    



void    AnalysisOfBetaUniverse::collectQFacesInWorlds(rg_dList<BetaFace*>& firstQFacesInWorlds) const
{
    //  for root world
    BetaVertex* virtualVertex = m_betaUniverse->findVirtualVertex();
    BetaCell*   firstCellOfVirtualVertex = virtualVertex->getFirstCell();
    rg_INT i = 0;
    for ( i=0; i<4; i++ ) {
        if ( firstCellOfVirtualVertex->getFace(i)->isVirtual() ) {
            firstQFacesInWorlds.add( virtualVertex->getFirstCell()->getFace(i) );
            break;
        }
    }

    rg_INT     numQTEdge = m_betaUniverse->getNumEdges();
    BetaEdge** qtEdge    = m_betaUniverse->getSortedEdges();

    for ( i=0; i<numQTEdge; i++ ) {
        if ( !qtEdge[i]->isGateToSmallWorlds() ) {
            continue;
        }

        rg_dList<BetaFace*>* facesInSmallWorld = qtEdge[i]->getSmallWorlds();
        facesInSmallWorld->reset4Loop();
        while ( facesInSmallWorld->setNext4Loop() ) {
            firstQFacesInWorlds.add( facesInSmallWorld->getEntity() );
        }
    }
}

