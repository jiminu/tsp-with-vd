// lusterFromQT.cpp: implementation of the SubQuasitriangulation class.
//
//////////////////////////////////////////////////////////////////////

#include "SubQuasitriangulation.h"
using namespace V::GeometryTier;



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SubQuasitriangulation::SubQuasitriangulation()
{    
    m_betaUniverse = rg_NULL;
}


SubQuasitriangulation::SubQuasitriangulation( BetaUniverse* betaUniverse )
{
    m_betaUniverse = betaUniverse;
}


SubQuasitriangulation::~SubQuasitriangulation()
{
}


void SubQuasitriangulation::makeClustersForSecondaryStructures(Molecule* aMolecule)
{   
    m_clusters.removeAll();

    rg_dList<Chain>* chains = aMolecule->getChains();

    chains->reset4Loop();
    while(chains->setNext4Loop()) {
        Chain* currChain = chains->getpEntity();

        if(currChain->getChainCode() != PROTEIN_CHAIN) {
            continue;
        }

        SecondaryStructure* currSecStruct = currChain->getSecondaryStructure();

        rg_dList<Helix>* helices = currSecStruct->getHelices();
        helices->reset4Loop();
        while(helices->setNext4Loop()) {
            Helix* currHelix = helices->getpEntity();
            ClusterFromQT* currCluster = m_clusters.add(ClusterFromQT(HELIX));

            rg_dList<BetaCell>* cells = m_betaUniverse->getCellList();
            cells->reset4Loop();
            while(cells->setNext4Loop()) {
                BetaCell* currCell = cells->getpEntity();

                if(currCell->isVirtual()) {
                    continue;
                }

                BetaVertex** vertices = currCell->getVertices();

                if(vertices[0]->isVirtual() || vertices[1]->isVirtual() || vertices[2]->isVirtual() || vertices[3]->isVirtual()) {
                    continue;
                }

                Atom* atoms[4] = {((Atom*)(vertices[0]->getProperty())),
                                  ((Atom*)(vertices[1]->getProperty())),
                                  ((Atom*)(vertices[2]->getProperty())),
                                  ((Atom*)(vertices[3]->getProperty())) };

                if(    currHelix->isResidueInHelix(atoms[0]->getResidue())
                    && currHelix->isResidueInHelix(atoms[1]->getResidue())
                    && currHelix->isResidueInHelix(atoms[2]->getResidue())
                    && currHelix->isResidueInHelix(atoms[3]->getResidue())) {

                        currCluster->addCellIntoHash(currCell);
                }
            }            
        }

        rg_dList<Sheet>* sheets = currSecStruct->getSheets();
        sheets->reset4Loop();
        while(sheets->setNext4Loop()) {
            Sheet* currSheet = sheets->getpEntity();
            Strand* strands = currSheet->getStrands();

            for(int i = 0; i < currSheet->getNumOfStrands(); i++) {                

                ClusterFromQT* currCluster = m_clusters.add(ClusterFromQT(BETASTRAND));

                rg_dList<BetaCell>* cells = m_betaUniverse->getCellList();
                cells->reset4Loop();
                while(cells->setNext4Loop()) {
                    BetaCell* currCell = cells->getpEntity();

                    if(currCell->isVirtual()) {
                        continue;
                    }

                    BetaVertex** vertices = currCell->getVertices();

                    if(vertices[0]->isVirtual() || vertices[1]->isVirtual() || vertices[2]->isVirtual() || vertices[3]->isVirtual()) {
                        continue;
                    }

                    Atom* atoms[4] = {((Atom*)(vertices[0]->getProperty())),
                                      ((Atom*)(vertices[1]->getProperty())),
                                      ((Atom*)(vertices[2]->getProperty())),
                                      ((Atom*)(vertices[3]->getProperty())) };

                    if(    strands[i].isResidueInStrand(atoms[0]->getResidue())
                        && strands[i].isResidueInStrand(atoms[1]->getResidue())
                        && strands[i].isResidueInStrand(atoms[2]->getResidue())
                        && strands[i].isResidueInStrand(atoms[3]->getResidue())) {

                            currCluster->addCellIntoHash(currCell);
                    }
                }       
            }
        }
    }
}



void SubQuasitriangulation::huntCellsInCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaCell* > >& betaCellList ) const
{
    m_clusters.reset4Loop();    
    while(m_clusters.setNext4Loop()) {
        rg_dList<BetaCell*>* currCellList = betaCellList.add(rg_dList<BetaCell*>());

        ClusterFromQT* currCluster = m_clusters.getpEntity();

        Hash4BetaCell* cells = currCluster->getBetaCells();

        Hash4BetaCell::iterator cell_i = cells->begin();

        for(; cell_i != cells->end(); cell_i++) {
            BetaCell* currCell = (*cell_i).second;

            if(currCell->getBoundingState(beta) != EXTERIOR_SIMPLEX) {
                currCellList->addWithoutSame(currCell);            
            }
        }        
    }
}


void SubQuasitriangulation::huntFacesInCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaFace* > >& betaFaceList ) const
{
    m_clusters.reset4Loop();    
    while(m_clusters.setNext4Loop()) {
        rg_dList<BetaFace*>* currFaceList = betaFaceList.add(rg_dList<BetaFace*>());

        ClusterFromQT* currCluster = m_clusters.getpEntity();
        
        Hash4BetaCell* cells = currCluster->getBetaCells();

        Hash4BetaCell::iterator cell_i = cells->begin();

        for(; cell_i != cells->end(); cell_i++) {
            BetaCell* currCell = (*cell_i).second;

            BetaFace** faces = currCell->getFaces();

            for(int i = 0; i < 4; i++) {
                if(faces[i]->getBoundingState(beta) != EXTERIOR_SIMPLEX) {
                    currFaceList->addWithoutSame(faces[i]);
                }
            }
        }        
    }
}

void SubQuasitriangulation::huntEdgesInCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaEdge* > >& betaEdgeList ) const
{
    m_clusters.reset4Loop();    
    while(m_clusters.setNext4Loop()) {
        rg_dList<BetaEdge*>* currEdgeList = betaEdgeList.add(rg_dList<BetaEdge*>());

        ClusterFromQT* currCluster = m_clusters.getpEntity();

        Hash4BetaCell* cells = currCluster->getBetaCells();

        Hash4BetaCell::iterator cell_i = cells->begin();

        for(; cell_i != cells->end(); cell_i++) {
            BetaCell* currCell = (*cell_i).second;

            BetaFace** faces = currCell->getFaces();            
            for(int i = 0; i < 4; i++) {
                rg_dList<BetaEdge*> boundingEdges;

                faces[i]->searchEdgesInIntraWorld(boundingEdges);

                boundingEdges.reset4Loop();
                while(boundingEdges.setNext4Loop()) {
                    BetaEdge* currEdge = boundingEdges.getEntity();

                    if(currEdge->getBoundingState(beta) != EXTERIOR_SIMPLEX) {
                        currEdgeList->addWithoutSame(currEdge);
                    }
                }
            }
        }        
    }
}

void SubQuasitriangulation::huntVerticesInCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaVertex* > >& betaVertexList ) const
{
    m_clusters.reset4Loop();    
    while(m_clusters.setNext4Loop()) {
        rg_dList<BetaVertex*>* currVertexList = betaVertexList.add(rg_dList<BetaVertex*>());

        ClusterFromQT* currCluster = m_clusters.getpEntity();

        Hash4BetaCell* cells = currCluster->getBetaCells();

        Hash4BetaCell::iterator cell_i = cells->begin();

        for(; cell_i != cells->end(); cell_i++) {
            BetaCell* currCell = (*cell_i).second;

            BetaFace** faces = currCell->getFaces();

            for(int i = 0; i < 4; i++) {
                BetaVertex* boundingVertex[3];
                faces[i]->searchVerticesInIntraWorld(boundingVertex);

                for(int j = 0; j < 3; j++) {
                    if(boundingVertex[j]->getBoundingState(beta) != EXTERIOR_SIMPLEX) {
                        currVertexList->addWithoutSame(boundingVertex[j]);
                    }
                }
            }
        }        
    }
}

void SubQuasitriangulation::huntFacesOnBoundaryOfCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaFace* > >& betaFaceList ) const
{
    rg_dList<ClusterFromQT> clusterList;
    huntCellsInCluster4GivenBetaValue( beta, clusterList );
    
    clusterList.reset4Loop();    
    while(clusterList.setNext4Loop()) {
        rg_dList<BetaFace*>* currFaceList = betaFaceList.add(rg_dList<BetaFace*>());

        ClusterFromQT* currCluster = clusterList.getpEntity();

        Hash4BetaCell* cells = currCluster->getBetaCells();

        Hash4BetaCell::iterator cell_i = cells->begin();

        for(; cell_i != cells->end(); cell_i++) {
            BetaCell* currCell = (*cell_i).second;

            BetaFace** faces = currCell->getFaces();

            for(int i = 0; i < 4; i++) {
                BetaCell* leftCell  = faces[i]->getLeftCell();
                BetaCell* rightCell = faces[i]->getRightCell();

                if(currCluster->isBetaCellInCluster(leftCell) == rg_TRUE && currCluster->isBetaCellInCluster(rightCell) == rg_TRUE) {
                    continue;
                }

                currFaceList->addWithoutSame(faces[i]);                
            }
        }        
    }
}

void SubQuasitriangulation::huntEdgesOnBoundaryOfCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaEdge* > >& betaEdgeList ) const
{
    rg_dList< rg_dList< BetaFace* > > betaFaceList;
    huntFacesOnBoundaryOfCluster4GivenBetaValue( beta, betaFaceList );

    betaFaceList.reset4Loop();
    while(betaFaceList.setNext4Loop()) {
        rg_dList<BetaEdge*>* currEdgeList = betaEdgeList.add(rg_dList<BetaEdge*>());

        rg_dList<BetaFace*>* currFaceList = betaFaceList.getpEntity();

        currFaceList->reset4Loop();
        while(currFaceList->setNext4Loop()) {
            BetaFace* currFace = currFaceList->getEntity();

            rg_dList<BetaEdge*> boundingEdges;

            currFace->searchEdgesInIntraWorld(boundingEdges);

            boundingEdges.reset4Loop();
            while(boundingEdges.setNext4Loop()) {
                BetaEdge* currEdge = boundingEdges.getEntity();
                currEdgeList->addWithoutSame(currEdge);            
            }
        }
    }
}

void SubQuasitriangulation::huntVerticesOnBoundaryOfCluster4GivenBetaValue( const rg_REAL& beta, rg_dList< rg_dList< BetaVertex* > >& betaVertexList ) const
{
    rg_dList< rg_dList< BetaFace* > > betaFaceList;
    huntFacesOnBoundaryOfCluster4GivenBetaValue( beta, betaFaceList );

    betaFaceList.reset4Loop();
    while(betaFaceList.setNext4Loop()) {
        rg_dList<BetaVertex*>* currVertexList = betaVertexList.add(rg_dList<BetaVertex*>());

        rg_dList<BetaFace*>* currFaceList = betaFaceList.getpEntity();

        currFaceList->reset4Loop();
        while(currFaceList->setNext4Loop()) {
            BetaFace* currFace = currFaceList->getEntity();

            BetaVertex* boundingVertex[3];
            currFace->searchVerticesInIntraWorld(boundingVertex);

            for(int j = 0; j < 3; j++) {                
                currVertexList->addWithoutSame(boundingVertex[j]);                
            }
        }
    }
}



void SubQuasitriangulation::huntCellsInCluster4GivenBetaValue( const rg_REAL& beta, rg_dList<ClusterFromQT>& clusterList ) const
{
    m_clusters.reset4Loop();    
    while(m_clusters.setNext4Loop()) {
        ClusterFromQT* currCluster = m_clusters.getpEntity();
        
        ClusterFromQT* currClusterForGivenBeta = clusterList.add(ClusterFromQT(UNK));
        currClusterForGivenBeta->setTyperOfCluster(currCluster->getTypeOfCluster());
        
        Hash4BetaCell* cells = currCluster->getBetaCells();

        Hash4BetaCell::iterator cell_i = cells->begin();

        for(; cell_i != cells->end(); cell_i++) {
            BetaCell* currCell = (*cell_i).second;

            if(currCell->getBoundingState(beta) != EXTERIOR_SIMPLEX) {
                currClusterForGivenBeta->addCellIntoHash(currCell);
            }
        }        
    }
}


ClusterType SubQuasitriangulation::getClusterType( const rg_INT& index ) const
{
    int indexCounter = 0;
    m_clusters.reset4Loop();
    while(m_clusters.setNext4Loop()) {
        if(indexCounter == index) {
            ClusterFromQT* currCluster = m_clusters.getpEntity();
            return currCluster->getTypeOfCluster();
        }
        else {
            indexCounter++;
        }
    }

    return UNK;
}