#include "AnalysisOfVoronoiDiagram.h"
#include "Bucket.h"
#include "BucketForSpheres.h"
using namespace V::GeometryTier;


AnalysisOfVoronoiDiagram::AnalysisOfVoronoiDiagram()
{
    m_VoronoiDiagram = rg_NULL;
}



AnalysisOfVoronoiDiagram::AnalysisOfVoronoiDiagram(rg_SphereSetVoronoiDiagram* diagram)
{
    m_VoronoiDiagram = diagram;
}



AnalysisOfVoronoiDiagram::AnalysisOfVoronoiDiagram(const AnalysisOfVoronoiDiagram& analysis)
{
    m_VoronoiDiagram = analysis.m_VoronoiDiagram;
}



AnalysisOfVoronoiDiagram::~AnalysisOfVoronoiDiagram()
{
}




AnalysisOfVoronoiDiagram& AnalysisOfVoronoiDiagram::operator =(const AnalysisOfVoronoiDiagram& analysis)
{
    if ( this == &analysis ) {
        return *this;
    }

    m_VoronoiDiagram = analysis.m_VoronoiDiagram;

    return *this;
}



void AnalysisOfVoronoiDiagram::reportCorrectnessOfVDGeometry(ofstream& fout)
{

    BucketForSpheres<void*> bucket2;
    rg_dList< pair<Sphere, void*> > balls;

    rg_dList< BallGenerator >* generators = m_VoronoiDiagram->getGlobalGeneratorList();
    generators->reset4Loop();
    while ( generators->setNext4Loop() ) {
        BallGenerator* currGen = generators->getpEntity();
        balls.add( make_pair( currGen->getBall(), currGen->getProperty() ) );
    }
    bucket2.construct( balls, m_VoronoiDiagram->getSizeOfBucketElement() );



}



void AnalysisOfVoronoiDiagram::reportCorrectnessOfVertexGeometry(ofstream& fout)
{
    BucketForSpheres<void*> bucket;
    rg_dList< pair<Sphere, void*> > balls;

    rg_dList< BallGenerator >* generators = m_VoronoiDiagram->getGlobalGeneratorList();
    generators->reset4Loop();
    while ( generators->setNext4Loop() ) {
        BallGenerator* currGen = generators->getpEntity();
        balls.add( make_pair( currGen->getBall(), currGen ) );
    }
    bucket.construct( balls, m_VoronoiDiagram->getSizeOfBucketElement() );



    int numBalls = 0;
    rg_dList< VDVertex >* voronoiVertices = m_VoronoiDiagram->getGlobalVerticeList();

    voronoiVertices->reset4Loop();
    while ( voronoiVertices->setNext4Loop() ) {
        VDVertex* currVtx = voronoiVertices->getpEntity();

        Sphere tangentSphere( currVtx->getPoint(), currVtx->getRadiusOfTangentSphere() );

        rg_dList<Sphere> sphereList;
        bucket.pickOutSpheresByMask( tangentSphere, sphereList );


        fout << currVtx->getID() << "\t";
        sphereList.reset4Loop();
        while ( sphereList.setNext4Loop() ) {
            Sphere currBall = sphereList.getEntity();

            rg_REAL dist = tangentSphere.distance( currBall );
            fout << dist << "\t";
            numBalls++;
        }
        fout << "\n";
    }
}    




void AnalysisOfVoronoiDiagram::reportCorrectnessOfEdgeGeometry(ofstream& fout)
{
}



