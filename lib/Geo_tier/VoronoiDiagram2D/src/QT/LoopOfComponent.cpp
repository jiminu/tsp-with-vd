//#include "StdAfx.h"
#include "LoopOfComponent.h"
#include "EdgeBU2D.h"
#include "VertexBU2D.h"
#include "ComponentBSFace.h"
using namespace V::GeometryTier;


LoopOfComponent::LoopOfComponent(void)
{
    m_isOuterLoop = false;
    m_length = 0.0;
}



LoopOfComponent::LoopOfComponent( const LoopOfComponent& loop )
{
    m_edgesOnLoop = loop.m_edgesOnLoop;
    m_isOuterLoop = loop.m_isOuterLoop;
    m_length      = loop.m_length;
}



LoopOfComponent::~LoopOfComponent(void)
{
}



double LoopOfComponent::computeLength()
{
    m_length = 0.0;

    m_edgesOnLoop.reset4Loop();
    while(m_edgesOnLoop.setNext4Loop()) {
        EdgeBU2D* currEdge = m_edgesOnLoop.getEntity();
        double    currLength = currEdge->getStartVertex()->getCoord().distance( currEdge->getEndVertex()->getCoord() );

        m_length += currLength;
    }

    return m_length;
}



double LoopOfComponent::computeSignedArea()
{
    double betaValue = m_component->getBetaValue();
    m_edgesOnLoop.reset4Loop();
    while ( m_edgesOnLoop.setNext4Loop() ) {
        EdgeBU2D* currEdge = m_edgesOnLoop.getEntity();
        rg_Point2D startPt, endPt;
        if ( currEdge->getLeftFace()->getBoundingState( betaValue ) == INTERIOR_SIMPLEX ) {
            startPt = currEdge->getStartVertex()->getCoord();
            endPt = currEdge->getEndVertex()->getCoord();
        }
        else {
            startPt = currEdge->getEndVertex()->getCoord();
            endPt = currEdge->getStartVertex()->getCoord();
        }
        
        double localSignedArea = 0.0;
        localSignedArea = ( startPt.getY() + endPt.getY() ) / 2.0 * ( startPt.getX() - endPt.getX() );
        m_signedArea += localSignedArea;
    }

    return m_signedArea;
}



void LoopOfComponent::findOrientedBoundingPoints( list<rg_Point2D>& boundingPoints )
{
    double betaValue = m_component->getBetaValue();

    m_edgesOnLoop.reset4Loop();
    while ( m_edgesOnLoop.setNext4Loop() ) {
        EdgeBU2D* edge = m_edgesOnLoop.getEntity();
        if ( edge->getLeftFace()->getBoundingState( betaValue ) == INTERIOR_SIMPLEX ) {
            boundingPoints.push_back( edge->getEndVertex()->getCoord() );
        }
        else {
            boundingPoints.push_back( edge->getStartVertex()->getCoord() );
        }
    }
}

void LoopOfComponent::findOrientedBoundingVertices(list<VertexBU2D*>& boundingVertices)
{
    double betaValue = m_component->getBetaValue();

    m_edgesOnLoop.reset4Loop();
    while (m_edgesOnLoop.setNext4Loop()) {
        EdgeBU2D* edge = m_edgesOnLoop.getEntity();
        if (edge->getLeftFace()->getBoundingState(betaValue) == INTERIOR_SIMPLEX) {
            boundingVertices.push_back(edge->getEndVertex());
        }
        else {
            boundingVertices.push_back(edge->getStartVertex());
        }
    }
}
