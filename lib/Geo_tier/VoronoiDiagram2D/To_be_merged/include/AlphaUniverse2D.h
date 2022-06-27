#ifndef _ALPHAUNIVERSE2D_H
#define _ALPHAUNIVERSE2D_H

#include "rg_dList.h"
#include "VertexBU2D.h"
#include "EdgeBU2D.h"
#include "FaceBU2D.h"
#include "Disc.h"

#include "Interval_REAL.h"

#include "DelaunayTriangulation2D_GEO.h"

namespace BULL2D {
namespace GeometryTier {


class AlphaUniverse2D
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
    AlphaUniverse2D();
    ~AlphaUniverse2D();

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


    void construct(DelaunayTriangulation2D_GEO& DT2D);
    void clean();
    void rearrangeSimplexIDs();


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
    void convert(DelaunayTriangulation2D_GEO& DT2D);
        void createAndInitializeFacesInBU2D( DelaunayTriangulation2D_GEO& DT2D,    map<FaceSCDS*, FaceBU2D*>&     mappingTableFromFaceInQTToFaceInBU );
        void createAndInitializeVerticesInBU2D( DelaunayTriangulation2D_GEO& DT2D, map<VertexSCDS*, VertexBU2D*>& mappingTableFromVtxInQTToVtxInBU );
    void computeBetaSpan();
    void generateSortedFiniteSimplexes();
        void generateSortedFiniteFaces();
        void generateSortedFiniteEdges();
        void generateSortedFiniteVertices();
        //friend int compareAlphaSpanOfFaces(const void* ts1, const void* ts2);
        //friend int compareAlphaSpanOfEdges(const void* ts1, const void* ts2);
        //friend int compareAlphaSpanOfVertices(const void* ts1, const void* ts2);
    void computeBetaInterval();
        void rearrangeFaceIDs();
        void rearrangeEdgeIDs();
        void rearrangeVertexIDs();
};


int compareAlphaSpanOfFaces(const void* ts1, const void* ts2);
int compareAlphaSpanOfEdges(const void* ts1, const void* ts2);
int compareAlphaSpanOfVertices(const void* ts1, const void* ts2);

} // namespace GeometryTier
} // namespace BULL2D

#endif

