#ifndef _BETAUNIVERSE2D_H
#define _BETAUNIVERSE2D_H

#include "rg_dList.h"
#include "VertexBU2D.h"
#include "EdgeBU2D.h"
#include "FaceBU2D.h"
#include "Disc.h"

#include "Interval_REAL.h"

#include "QuasiTriangulation2D.h"

#include <vector>
using namespace std;


namespace V {
namespace GeometryTier {


class BetaUniverse2D
{
private:
    //  extended quasi-simplicial complex data structure
    rg_dList<VertexBU2D> m_vertices;
    rg_dList<EdgeBU2D>   m_edges;
    rg_dList<FaceBU2D>   m_faces;

    rg_dList<Disc>       m_discs;

    //  if beta <= m_betaInterval.m_lowerValue, its beta-shape is equal to the set of center points of discs.
    //  if beta >= m_betaInterval.m_upperValue, its beta-shape is equal to squeezed-Hull.
    Interval_REAL        m_betaInterval;

	// sorted simplexes
	VertexBU2D**         m_sortedFiniteVertices;
	EdgeBU2D**           m_sortedFiniteEdges;
	FaceBU2D**           m_sortedFiniteFaces;

	rg_INT               m_numFiniteVertices;
	rg_INT               m_numFiniteEdges;
	rg_INT               m_numFiniteFaces;

public:
    BetaUniverse2D();
    ~BetaUniverse2D();

    inline rg_dList<VertexBU2D>&   getVertices() {return m_vertices;}
    inline rg_dList<EdgeBU2D>&     getEdges()    {return m_edges;}
    inline rg_dList<FaceBU2D>&     getFaces()    {return m_faces;}
    inline rg_dList<Disc>&         getDiscs()    {return m_discs;}

	inline Interval_REAL&          getBetaInterval() {return m_betaInterval;}

    inline rg_INT getNumVertices() const         {return m_numFiniteVertices;}
    inline rg_INT getNumEdges() const            {return m_numFiniteEdges;}
    inline rg_INT getNumFaces() const            {return m_numFiniteFaces;}
    inline rg_INT getNumDiscs() const            {return m_discs.getSize();}

    inline VertexBU2D** getSortedVertices()      {return m_sortedFiniteVertices;}
    inline EdgeBU2D**   getSortedEdges()         {return m_sortedFiniteEdges;}
    inline FaceBU2D**   getSortedFaces()         {return m_sortedFiniteFaces;}


    void construct(QuasiTriangulation2D& QT2D);
    void clean();
    void rearrangeSimplexIDs();

    // Mokwon added 2017.10.24 for image processing.
    void        constructForImageProcessing(const vector< vector< vector< int > > >& colorDataSet, const int& lengthOfXDir, const int& lengthOfYDir, const bool& darkerPixelBecomesBiggerDisk );
    rg_FLOAT    calculateRadiusOfPixel( const int& r, const int& g, const int& b, const bool& darkerPixelBecomesBiggerDisk );
    bool        diagonal_is_positive( const Disc* const disk_leftBottom, const Disc* const disk_rightBottom, const Disc* const disk_leftTop, const Disc* const disk_rightTop );
    void        updateBetaSpan( const int& lengthOfX, const int& lengthOfY );
    rg_FLOAT    calculateRadiusOfTangentCircle_quadraticFormula(const double& baseRadius, const double& prevRadius_CCW, const double& nextRadius_CCW);
    void        computeBetaSpanForImageProcessing( const int& lengthOfX, const int& lengthOfY );
    // Mokwon added Mar 9, 2021 for Hyunwoo.
    void        constructForImageProcessing( const vector< vector< int > >& grayDataSet, const bool& darkerPixelBecomesBiggerDisk );
    // Joonghyun added July 23, 2020 (tentative)
    void        constructForImageProcessing(const vector< vector< unsigned char > >& colorDataSet, const int& lengthOfXDir, const int& lengthOfYDir, const bool& darkerPixelBecomesBiggerDisk);
    rg_FLOAT    calculateRadiusOfPixel(const unsigned char& intensity, const bool& darkerPixelBecomesBiggerDisk);
    ////////////////////////////////////////////////

	// Query functions
	void  huntFacesInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const;
	void  huntEdgesInBetaComplex(    const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const;
	void  huntVerticesInBetaComplex( const rg_REAL& beta, rg_INT& startIndex, rg_INT& endIndex ) const;

	void  huntFacesInBetaComplex(    const rg_REAL& beta, rg_dList<FaceBU2D*>&   betaFaces ) const;
	void  huntEdgesInBetaComplex(    const rg_REAL& beta, rg_dList<EdgeBU2D*>&   betaEdges ) const;
	void  huntVerticesInBetaComplex( const rg_REAL& beta, rg_dList<VertexBU2D*>& betaVertices ) const;

	void  huntEdgesInBetaShape(    const rg_REAL& beta, rg_dList<EdgeBU2D*>&   betaEdges ) const;
	void  huntVerticesInBetaShape( const rg_REAL& beta, rg_dList<VertexBU2D*>& betaVertices ) const;


    void    reportBetaUniverse2D(const string& reportFilename);


private:
    void convert(QuasiTriangulation2D& QT2D);
        void createAndInitializeFacesInBU2D( QuasiTriangulation2D& QT2D,    map<FaceSCDS*, FaceBU2D*>&     mappingTableFromFaceInQTToFaceInBU );
        void createAndInitializeVerticesInBU2D( QuasiTriangulation2D& QT2D, map<VertexSCDS*, VertexBU2D*>& mappingTableFromVtxInQTToVtxInBU );
    void computeBetaSpan();
    void generateSortedFiniteSimplexes();
        void generateSortedFiniteFaces();
        void generateSortedFiniteEdges();
        void generateSortedFiniteVertices();
        friend int compareBetaSpanOfFaces(const void* ts1, const void* ts2);
        friend int compareBetaSpanOfEdges(const void* ts1, const void* ts2);
        friend int compareBetaSpanOfVertices(const void* ts1, const void* ts2);
    void computeBetaInterval();
        void rearrangeFaceIDs();
        void rearrangeEdgeIDs();
        void rearrangeVertexIDs();
};


} // GeometryTier
} // V

#endif

