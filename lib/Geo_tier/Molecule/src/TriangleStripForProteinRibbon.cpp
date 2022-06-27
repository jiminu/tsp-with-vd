#include "TriangleStripForProteinRibbon.h"
#include "rg_Atom.h"
#include "Chain.h"
#include "rg_BSplineCurve3D.h"
using namespace V::GeometryTier;

TriangleStripForProteinRibbon::TriangleStripForProteinRibbon()
{
    m_numVertices = 0;
    m_vertex[0]   = rg_NULL;
    m_vertex[1]   = rg_NULL;
    m_normal[0]   = rg_NULL;
    m_normal[1]   = rg_NULL;
}



TriangleStripForProteinRibbon::TriangleStripForProteinRibbon(const TriangleStripForProteinRibbon& tStrip)
{
    m_numVertices = tStrip.m_numVertices;

    for ( int i=0; i<2; i++)  {
        m_vertex[i] = new rg_Point3D[m_numVertices];
        m_normal[i] = new rg_Point3D[m_numVertices];

        for (int j=0; j<m_numVertices; j++)  {
            m_vertex[i][j] = tStrip.m_vertex[i][j];
            m_normal[i][j] = tStrip.m_normal[i][j];
        }
    }
}



TriangleStripForProteinRibbon::~TriangleStripForProteinRibbon()
{
    for ( int i=0; i<2; i++)  {
        if ( m_vertex[i] != rg_NULL )
            delete [] m_vertex[i];
        if ( m_normal[i] != rg_NULL )
            delete [] m_normal[i];
    }
}


    
void TriangleStripForProteinRibbon::clean()
{
    for ( rg_INT i=0; i<2; i++)  {
        if ( m_vertex[i] != rg_NULL ) {
            delete [] m_vertex[i];
        }

        if ( m_normal[i] != rg_NULL ) {
            delete [] m_normal[i];
        }
    }
    m_numVertices = 0;
}


rg_INT      TriangleStripForProteinRibbon::getNumVertices() const
{
    return m_numVertices;
}



rg_Point3D** TriangleStripForProteinRibbon::getVertices()
{
    return m_vertex;
}



rg_Point3D** TriangleStripForProteinRibbon::getNormals()
{
    return m_normal;
}




void TriangleStripForProteinRibbon::makeProteinRibbonByCarsonBugg86(Chain* currChain)
{
    if ( currChain->getChainCode() != PROTEIN_CHAIN )  {
        return;
    }

    clean();



    //  get residues in current chain.
    rg_dList<Residue*> residuesInCurrChain;
	currChain->getResiduesSortedBySequenceNumber( residuesInCurrChain );

    rg_INT numResidues = residuesInCurrChain.getSize();
    rg_INT residueID   = 0;
    Residue** residue  = new Residue*[numResidues];

    residuesInCurrChain.reset4Loop();
    while ( residuesInCurrChain.setNext4Loop() ) {
        Residue* currRes = residuesInCurrChain.getEntity();

        Atom* alphaCarbon    = currRes->getAlphaCarbonOfAminoResidue();
        Atom* carbonylOxygen = currRes->getOxygenInCarboxylGroupOfAminoResidue();
		if ( alphaCarbon != rg_NULL && carbonylOxygen != rg_NULL ) {
			residue[residueID] = currRes;
			residueID++;
		}
    }
    rg_INT       numResiduesInCurrChain = residueID;



    //  define two guide-curves
    rg_INT       i = 0;
    const rg_INT numGuideCurve = 2;
    rg_INT       numCtrlPtsOfGuideCurve = numResiduesInCurrChain+1; 
    rg_Point3D*  ctrlPt[numGuideCurve];
    for ( i=0; i<numGuideCurve; i++) {
        ctrlPt[i]    = new rg_Point3D[ numCtrlPtsOfGuideCurve ];
    }


    //  set control points of two guide-curves
    SecondaryStructure* secondaryStruct = currChain->getSecondaryStructure();

    rg_Point3D prevVectorForWidth;
    rg_INT i_res = 0;
    for ( i_res=0; i_res<numResiduesInCurrChain-1; i_res++)  {
        Atom* alphaCarbon     = residue[i_res]->getAlphaCarbonOfAminoResidue();
        Atom* nextAlphaCarbon = residue[i_res+1]->getAlphaCarbonOfAminoResidue();
        Atom* carbonylOxygen  = residue[i_res]->getOxygenInCarboxylGroupOfAminoResidue();

        rg_Point3D currCarbon = alphaCarbon->getpAtomBall()->getCenter();
        rg_Point3D nextCarbon = nextAlphaCarbon->getpAtomBall()->getCenter();
        rg_Point3D currOxygen = carbonylOxygen->getpAtomBall()->getCenter();

        rg_Point3D vecA = nextCarbon - currCarbon;
        rg_Point3D vecB = currOxygen - currCarbon;
        rg_Point3D normalOfPeptidePlane = vecA.crossProduct(vecB);
        rg_Point3D vectorForWidth       = normalOfPeptidePlane.crossProduct(vecA);

        normalOfPeptidePlane.normalize();
        vectorForWidth.normalize();

        rg_Point3D pt = (currCarbon + nextCarbon)/2.0;
        //  produce a reasonable helix diameter.
        if ( secondaryStruct->isResidueInHelix(residue[i_res]) ) {
            pt += 1.5*normalOfPeptidePlane;
        }

        //  determine ribbon width 
        rg_REAL ribbonWidth = 0.0;
        if (    secondaryStruct->isResidueInHelix(residue[i_res]) 
             || secondaryStruct->isResidueInSheet(residue[i_res]) ) {
            ribbonWidth = 1.5;
        }
        else {
            ribbonWidth = 0.5;
        }

        //  determine control points of curve.
        if ( i_res == 0 ) {
            prevVectorForWidth = vectorForWidth;
        }
        rg_REAL angle = vectorForWidth.angle( prevVectorForWidth );
        if ( angle < ANGLE_90_DEGREE ) {            
            ctrlPt[0][i_res+1] = pt + (ribbonWidth*vectorForWidth);
            ctrlPt[1][i_res+1] = pt - (ribbonWidth*vectorForWidth);
            prevVectorForWidth = vectorForWidth;
        }
        else {
            ctrlPt[0][i_res+1] = pt - (ribbonWidth*vectorForWidth);
            ctrlPt[1][i_res+1] = pt + (ribbonWidth*vectorForWidth);
            prevVectorForWidth = -vectorForWidth;
        }
    }
    
    rg_Point3D firstCtrlPt = residue[0]->getAlphaCarbonOfAminoResidue()->getpAtomBall()->getCenter();
    rg_Point3D lastCtrlPt  = residue[numResiduesInCurrChain-1]->getAlphaCarbonOfAminoResidue()->getpAtomBall()->getCenter();
    ctrlPt[1][0] = ctrlPt[0][0] = firstCtrlPt;
    ctrlPt[1][numCtrlPtsOfGuideCurve-1] = ctrlPt[0][numCtrlPtsOfGuideCurve-1] = lastCtrlPt;



    //  make guide-curves and get points on guide-curves
    rg_INT numLinePerCurveSegment = 10;
    m_numVertices = numCtrlPtsOfGuideCurve*numLinePerCurveSegment;

    rg_BSplineCurve3D guideCurve[numGuideCurve];
    rg_BSplineCurve3D derivativeOfGuideCurve[numGuideCurve];

    for ( i=0; i<numGuideCurve; i++)  {
        guideCurve[i].setOrder( 4 );
        guideCurve[i].setCtrlPts( numCtrlPtsOfGuideCurve, ctrlPt[i] );
        delete [] ctrlPt[i];

        derivativeOfGuideCurve[i] = guideCurve[i].makeDerivative();

        m_vertex[i] = new rg_Point3D[ m_numVertices ];
        m_normal[i] = new rg_Point3D[ m_numVertices ];

    }

    rg_REAL maxParamINT = (rg_REAL)(m_numVertices-1);
    for ( rg_INT i_u = 0; i_u<m_numVertices; i_u++ ) {
        rg_REAL u = i_u/maxParamINT;

        m_vertex[0][i_u] = guideCurve[0].evaluatePt( u );
        m_vertex[1][i_u] = guideCurve[1].evaluatePt( u );

        rg_Point3D tangent[3];
        tangent[0] = derivativeOfGuideCurve[0].evaluatePt( u );
        tangent[1] = derivativeOfGuideCurve[1].evaluatePt( u );
        tangent[2] = m_vertex[1][i_u] - m_vertex[0][i_u];

        m_normal[0][i_u] = tangent[2].crossProduct(tangent[0]);
        m_normal[0][i_u].normalize();
        m_normal[1][i_u] = tangent[2].crossProduct(tangent[1]);
        m_normal[1][i_u].normalize();

    }
}



TriangleStripForProteinRibbon& TriangleStripForProteinRibbon::operator =(const TriangleStripForProteinRibbon& tStrip)
{
    if ( this == &tStrip )
        return *this;

    m_numVertices = tStrip.m_numVertices;

    for ( int i=0; i<2; i++)  {
        m_vertex[i] = new rg_Point3D[m_numVertices];
        m_normal[i] = new rg_Point3D[m_numVertices];

        for (int j=0; j<m_numVertices; j++)  {
            m_vertex[i][j] = tStrip.m_vertex[i][j];
            m_normal[i][j] = tStrip.m_normal[i][j];
        }
    }

    return *this;
}



