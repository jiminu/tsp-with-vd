#ifndef _RG_WINGEDEDGEDATASTRUCTURE_H_
#define _RG_WINGEDEDGEDATASTRUCTURE_H_

#include "rg_Const.h"
#include "rg_dList.h"

template<class VertexGeom, class EdgeGeom, class FaceGeom>
class WVertex;

template<class VertexGeom, class EdgeGeom, class FaceGeom>
class WFace;

template<class VertexGeom, class EdgeGeom, class FaceGeom>
class WEdge
{
private:
    WVertex<VertexGeom,EdgeGeom,FaceGeom>* m_startVertex;
    WVertex<VertexGeom,EdgeGeom,FaceGeom>* m_endVertex;
    
    WFace<VertexGeom,EdgeGeom,FaceGeom>*   m_rightFace;
    WFace<VertexGeom,EdgeGeom,FaceGeom>*   m_leftFace;

    WEdge<VertexGeom,EdgeGeom,FaceGeom>*   m_rightHand;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>*   m_leftHand;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>*   m_rightLeg;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>*   m_leftLeg;

	rg_INT	 m_ID;
	EdgeGeom m_geometry;

    rg_BOOL  m_visited;

public:


    //constructors
    WEdge();
    WEdge(const rg_INT& ID);
    WEdge(const rg_INT& tID,
          WVertex<VertexGeom,EdgeGeom,FaceGeom>* tStartVertex,
          WVertex<VertexGeom,EdgeGeom,FaceGeom>* tEndVertex);
    WEdge(WVertex<VertexGeom,EdgeGeom,FaceGeom>*      tStartVertex,
               WVertex<VertexGeom,EdgeGeom,FaceGeom>* tEndVertex,
               WFace<VertexGeom,EdgeGeom,FaceGeom>*   tRightFace,
               WFace<VertexGeom,EdgeGeom,FaceGeom>*   tLeftFace,
               WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightHand,
               WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftHand,
               WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightLeg,
               WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftLeg,
			   const rg_INT&	tID);
    WEdge(const WEdge<VertexGeom,EdgeGeom,FaceGeom>& tempEdge);

    ~WEdge();

    //get functions
    inline WVertex<VertexGeom,EdgeGeom,FaceGeom>* getStartVertex() const { return m_startVertex; }
    inline WVertex<VertexGeom,EdgeGeom,FaceGeom>* getEndVertex() const   { return m_endVertex; }

    inline WFace<VertexGeom,EdgeGeom,FaceGeom>*   getRightFace() const   { return m_rightFace; }
    inline WFace<VertexGeom,EdgeGeom,FaceGeom>*   getLeftFace() const    { return m_leftFace; }

    inline WEdge<VertexGeom,EdgeGeom,FaceGeom>*   getRightHand() const   { return m_rightHand; }
    inline WEdge<VertexGeom,EdgeGeom,FaceGeom>*   getLeftHand() const    { return m_leftHand; }
    inline WEdge<VertexGeom,EdgeGeom,FaceGeom>*   getRightLeg() const    { return m_rightLeg; }
    inline WEdge<VertexGeom,EdgeGeom,FaceGeom>*   getLeftLeg() const     { return m_leftLeg; }

	rg_INT			getID() const;
	EdgeGeom		getGeometry() const;
	EdgeGeom*		getpGeometry();

    inline rg_BOOL isVisited() const { return m_visited; }
    void    isVisited(const rg_BOOL& visited);

    //set functions
    void setStartVertex(WVertex<VertexGeom,EdgeGeom,FaceGeom>* tempVertex);
    void setEndVertex(WVertex<VertexGeom,EdgeGeom,FaceGeom>* tempVertex);

    void setRightFace(WFace<VertexGeom,EdgeGeom,FaceGeom>* tempFace);
    void setLeftFace(WFace<VertexGeom,EdgeGeom,FaceGeom>* tempFace);

    void setRightHand(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge);
    void setLeftHand(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge);
    void setRightLeg(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge);
    void setLeftLeg(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge);

	void setID(const rg_INT& tID);

	void setEdge(WFace<VertexGeom,EdgeGeom,FaceGeom>*   tRightFace,
                 WFace<VertexGeom,EdgeGeom,FaceGeom>*   tLeftFace,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightHand,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftHand,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightLeg,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftLeg);
	void setEdge(WVertex<VertexGeom,EdgeGeom,FaceGeom>* tStartVertex,
                 WVertex<VertexGeom,EdgeGeom,FaceGeom>* tEndVertex,
                 WFace<VertexGeom,EdgeGeom,FaceGeom>*   tRightFace,
                 WFace<VertexGeom,EdgeGeom,FaceGeom>*   tLeftFace,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightHand,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftHand,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightLeg,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftLeg );
	void setEdge(WVertex<VertexGeom,EdgeGeom,FaceGeom>* tStartVertex,
                 WVertex<VertexGeom,EdgeGeom,FaceGeom>* tEndVertex,
                 WFace<VertexGeom,EdgeGeom,FaceGeom>*   tRightFace,
                 WFace<VertexGeom,EdgeGeom,FaceGeom>*   tLeftFace,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightHand,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftHand,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightLeg,
                 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftLeg,
				 const rg_INT&	  tID);

	void setGeometry(const EdgeGeom& tGeometry);

    bool    isThere( WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex );
    inline  rg_BOOL isIncidentTo(WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex) const {
                if ( m_startVertex == vertex || m_endVertex == vertex ) return true;
                else return false; 
            }

    rg_BOOL isConnected(WEdge<VertexGeom,EdgeGeom,FaceGeom>* otherEdge);

    //operator
    WEdge<VertexGeom,EdgeGeom,FaceGeom>& operator=(const WEdge<VertexGeom,EdgeGeom,FaceGeom>& tempEdge);
    bool operator==(const WEdge<VertexGeom,EdgeGeom,FaceGeom>& tempEdge) const;

    WVertex<VertexGeom,EdgeGeom,FaceGeom>* getVertexOfLeftHand() const;
    WVertex<VertexGeom,EdgeGeom,FaceGeom>* getVertexOfRightHand() const;

    void connectRightHand(WEdge<VertexGeom,EdgeGeom,FaceGeom>* rightHand);
    void connectLeftHand( WEdge<VertexGeom,EdgeGeom,FaceGeom>* leftHand);
    void connectRightLeg( WEdge<VertexGeom,EdgeGeom,FaceGeom>* rightLeg);
    void connectLeftLeg(  WEdge<VertexGeom,EdgeGeom,FaceGeom>* leftLeg);

    void connectRightFace( WFace<VertexGeom,EdgeGeom,FaceGeom>* rightFace);
    void connectLeftFace(  WFace<VertexGeom,EdgeGeom,FaceGeom>* leftFace);


    void reverse();


    void getEdgesInStar(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList);
    void getFacesInStar(rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*>& faceList);

};


template<class VertexGeom, class EdgeGeom, class FaceGeom>
class WFace
{
private:
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* m_firstEdge;

	rg_INT		m_ID;
    FaceGeom    m_geometry;

    rg_BOOL     m_visited;
    
public:
    //constructors
    WFace();
    WFace(const  rg_INT&		tID);
    WFace(const  rg_INT&		tID,
          WEdge<VertexGeom,EdgeGeom,FaceGeom>* tFirstEdge);
    WFace(       WEdge<VertexGeom,EdgeGeom,FaceGeom>* tFirstEdge,
          const  FaceGeom&   tGeometry,
	      const  rg_INT&		tID);
    WFace(const  WFace<VertexGeom,EdgeGeom,FaceGeom>&       tempFace);
    ~WFace();

    //get functions
    inline WEdge<VertexGeom,EdgeGeom,FaceGeom>* getFirstEdge() const { return m_firstEdge; }
    FaceGeom		getGeometry() const;
    FaceGeom*		getpGeometry();
	rg_INT			getID() const;

    inline rg_BOOL isVisited() const { return m_visited; }
    void    isVisited(const rg_BOOL& visited);

    rg_BOOL isThere(WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex);
    rg_BOOL isIncidentTo(WEdge<VertexGeom,EdgeGeom,FaceGeom>* edge) const;

    //set functions
    void setFirstEdge(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge);
    void setGeometry(const FaceGeom& tGeometry);
	void setID(const rg_INT& tID);

    WFace& operator=(const WFace<VertexGeom,EdgeGeom,FaceGeom>& tempFace);

	rg_FLAG getIncidentVertexList(rg_dList<WVertex<VertexGeom,EdgeGeom,FaceGeom>*>& vertexList);
	rg_FLAG getIncidentEdgeList(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList);
	rg_FLAG getIncidentFaceList(rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*>& faceList);

	rg_FLAG getBoundingVertexList(rg_dList<WVertex<VertexGeom,EdgeGeom,FaceGeom>*>& boundingVertices);
    rg_FLAG getBoundingEdgeList(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& boundingEdges);


};


template<class VertexGeom,class EdgeGeom,class FaceGeom>
class WVertex
{
private:
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* m_firstEdge;
	rg_INT			m_ID;
	VertexGeom		m_geometry;
	
    rg_BOOL         m_visited;

public:
    //consturctors
    WVertex();
    WVertex(const rg_INT& tID);
    WVertex(const rg_INT& tID, const VertexGeom& tGeometry);
    WVertex(const WVertex<VertexGeom,EdgeGeom,FaceGeom>& tempVertex);

    ~WVertex();

    //get functions
    inline WEdge<VertexGeom,EdgeGeom,FaceGeom>* getFirstEdge()  const { return m_firstEdge; }
    
    inline VertexGeom	getGeometry() const { return m_geometry; }
	VertexGeom*	getpGeometry();	
	rg_INT		getID()			const;

    inline rg_BOOL isVisited() const { return m_visited; }
    void    isVisited(const rg_BOOL& visited);


    //set functions
    void setFirstEdge(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge);
    void setGeometry(const VertexGeom& tGeometry);
	void setID(const rg_INT& tID);


    //operator overloading
    WVertex<VertexGeom,EdgeGeom,FaceGeom>& operator=(const WVertex<VertexGeom,EdgeGeom,FaceGeom>& tempVertex);
    bool operator==(const WVertex<VertexGeom,EdgeGeom,FaceGeom>& tempVertex) const;

	//Query
	rg_FLAG getNeighborWVertexList(rg_dList<WVertex<VertexGeom,EdgeGeom,FaceGeom>*>& vertexList);
	rg_FLAG getIncidentWEdgeList(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList);
	rg_FLAG getIncidentWFaceList(rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*>& faceList);

    void    getEdgesInStar(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList);
    void    getEdgesInShell(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList);

    WEdge<VertexGeom,EdgeGeom,FaceGeom>* findConnectingEdge(WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex);
    rg_INT                               findSharingFace(WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex,
                                                         rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*>& faceList);

};


	




template<class VertexGeom, class EdgeGeom, class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>::WEdge()
{
    m_startVertex = (WVertex<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_endVertex   = (WVertex<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;

    m_rightFace = (WFace<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftFace  = (WFace<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;

    m_rightHand = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftHand  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_rightLeg  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftLeg  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;

	m_ID = -1;
    m_visited  = rg_FALSE;

}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>::WEdge(const rg_INT& ID)
{
	m_ID = ID;
    m_startVertex = (WVertex<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_endVertex   = (WVertex<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;

    m_rightFace = (WFace<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftFace  = (WFace<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;

    m_rightHand = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftHand  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_rightLeg  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftLeg  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;

    m_visited  = rg_FALSE;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>::WEdge(const rg_INT& tID,
                                           WVertex<VertexGeom,EdgeGeom,FaceGeom>* tStartVertex,
                                           WVertex<VertexGeom,EdgeGeom,FaceGeom>* tEndVertex)
{
	m_ID          = tID;
    m_startVertex = tStartVertex;
    m_endVertex   = tEndVertex;

    m_rightFace = (WFace<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftFace  = (WFace<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;

    m_rightHand = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftHand  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_rightLeg  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_leftLeg  = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;

    m_visited  = rg_FALSE;

}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>::WEdge(WVertex<VertexGeom,EdgeGeom,FaceGeom>*      tStartVertex,
                                           WVertex<VertexGeom,EdgeGeom,FaceGeom>*      tEndVertex,
                                           WFace<VertexGeom,EdgeGeom,FaceGeom>*        tRightFace,
                                           WFace<VertexGeom,EdgeGeom,FaceGeom>*        tLeftFace,
                                           WEdge<VertexGeom,EdgeGeom,FaceGeom>*  tRightHand,
                                           WEdge<VertexGeom,EdgeGeom,FaceGeom>*  tLeftHand,
                                           WEdge<VertexGeom,EdgeGeom,FaceGeom>*  tRightLeg,
                                           WEdge<VertexGeom,EdgeGeom,FaceGeom>*  tLeftLeg,
					                       const rg_INT&	tID)
{
    m_startVertex = tStartVertex;
    m_endVertex   = tEndVertex;
    m_rightFace   = tRightFace;
    m_leftFace    = tLeftFace;
    m_rightHand   = tRightHand;
    m_leftHand    = tLeftHand;
    m_rightLeg    = tRightLeg;
    m_leftLeg     = tLeftLeg;
	m_ID          = tID;
    m_visited  = rg_FALSE;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>::WEdge(const WEdge<VertexGeom,EdgeGeom,FaceGeom>& tempEdge)
{
    m_startVertex = tempEdge.m_startVertex;
    m_endVertex   = tempEdge.m_endVertex;

    m_rightFace   = tempEdge.m_rightFace;
    m_leftFace    = tempEdge.m_leftFace;

    m_rightHand   = tempEdge.m_rightHand;
    m_leftHand    = tempEdge.m_leftHand;
    m_rightLeg    = tempEdge.m_rightLeg;
    m_leftLeg     = tempEdge.m_leftLeg;

    m_geometry    = tempEdge.m_geometry;
	m_ID          = tempEdge.m_ID;
    m_visited     = tempEdge.m_visited;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>::~WEdge()
{
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
rg_INT WEdge<VertexGeom,EdgeGeom,FaceGeom>::getID() const
{
	return m_ID;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
EdgeGeom WEdge<VertexGeom,EdgeGeom,FaceGeom>::getGeometry() const
{
	return m_geometry;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
EdgeGeom* WEdge<VertexGeom,EdgeGeom,FaceGeom>::getpGeometry() 
{
	return &m_geometry;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::isVisited(const rg_BOOL& visited)
{
    m_visited = visited;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setStartVertex(WVertex<VertexGeom,EdgeGeom,FaceGeom>* tempVertex)
{
    m_startVertex = tempVertex;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setEndVertex(WVertex<VertexGeom,EdgeGeom,FaceGeom>* tempVertex)
{
    m_endVertex = tempVertex;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setRightFace(WFace<VertexGeom,EdgeGeom,FaceGeom>* tempFace)
{
    m_rightFace=tempFace;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setLeftFace(WFace<VertexGeom,EdgeGeom,FaceGeom>* tempFace)
{
    m_leftFace=tempFace;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setRightHand(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge)
{
    m_rightHand=tempEdge;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setLeftHand(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge)
{
    m_leftHand=tempEdge;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setRightLeg(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge)
{
    m_rightLeg=tempEdge;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setLeftLeg(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge)
{
    m_leftLeg=tempEdge;
}


template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setEdge(WFace<VertexGeom,EdgeGeom,FaceGeom>*   tRightFace,
                                                  WFace<VertexGeom,EdgeGeom,FaceGeom>*   tLeftFace,
                                                  WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightHand,
                                                  WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftHand,
                                                  WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightLeg,
                                                  WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftLeg)
{
    m_rightFace   = tRightFace;
    m_leftFace    = tLeftFace;
    m_rightHand   = tRightHand;
    m_leftHand    = tLeftHand;
    m_rightLeg    = tRightLeg;
    m_leftLeg     = tLeftLeg;
}




template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setEdge(
                         WVertex<VertexGeom,EdgeGeom,FaceGeom>* tStartVertex,
						 WVertex<VertexGeom,EdgeGeom,FaceGeom>* tEndVertex,
					     WFace<VertexGeom,EdgeGeom,FaceGeom>*   tRightFace,
						 WFace<VertexGeom,EdgeGeom,FaceGeom>*   tLeftFace,
						 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightHand,
						 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftHand,
						 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tRightLeg,
						 WEdge<VertexGeom,EdgeGeom,FaceGeom>*   tLeftLeg)
{
    m_startVertex = tStartVertex;
    m_endVertex   = tEndVertex;
    m_rightFace   = tRightFace;
    m_leftFace    = tLeftFace;
    m_rightHand   = tRightHand;
    m_leftHand    = tLeftHand;
    m_rightLeg    = tRightLeg;
    m_leftLeg     = tLeftLeg;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setEdge(WVertex<VertexGeom,EdgeGeom,FaceGeom>*      tStartVertex,
						 WVertex<VertexGeom,EdgeGeom,FaceGeom>*      tEndVertex,
					     WFace<VertexGeom,EdgeGeom,FaceGeom>*        tRightFace,
						 WFace<VertexGeom,EdgeGeom,FaceGeom>*        tLeftFace,
						 WEdge<VertexGeom,EdgeGeom,FaceGeom>*  tRightHand,
						 WEdge<VertexGeom,EdgeGeom,FaceGeom>*  tLeftHand,
						 WEdge<VertexGeom,EdgeGeom,FaceGeom>*  tRightLeg,
						 WEdge<VertexGeom,EdgeGeom,FaceGeom>*  tLeftLeg,
						 const rg_INT&	  tID)
{
    m_startVertex = tStartVertex;
    m_endVertex   = tEndVertex;
    m_rightFace   = tRightFace;
    m_leftFace    = tLeftFace;
    m_rightHand   = tRightHand;
    m_leftHand    = tLeftHand;
    m_rightLeg    = tRightLeg;
    m_leftLeg     = tLeftLeg;
	m_ID = tID;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setID(const rg_INT& tID)
{
	m_ID = tID;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::setGeometry(const EdgeGeom& tGeometry)
{
	m_geometry = tGeometry;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
bool WEdge<VertexGeom,EdgeGeom,FaceGeom>::isThere( WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex )
{
    if ( m_startVertex == vertex || m_endVertex == vertex )
        return true;
    else
        return false;
}




template<class VertexGeom, class EdgeGeom, class FaceGeom>
rg_BOOL WEdge<VertexGeom,EdgeGeom,FaceGeom>::isConnected(WEdge<VertexGeom,EdgeGeom,FaceGeom>* otherEdge)
{
    if ( m_startVertex == otherEdge->m_startVertex ) {
        return rg_TRUE;
    }
    else if ( m_startVertex == otherEdge->m_endVertex ) {
        return rg_TRUE;
    }
    else if ( m_endVertex == otherEdge->m_startVertex ) {
        return rg_TRUE;
    }
    else if ( m_endVertex == otherEdge->m_endVertex ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>& WEdge<VertexGeom,EdgeGeom,FaceGeom>::operator=(const WEdge<VertexGeom,EdgeGeom,FaceGeom>& tempEdge)
{
    if(this==&tempEdge)
    {
        return *this;
    }

    m_startVertex = tempEdge.m_startVertex;
    m_endVertex   = tempEdge.m_endVertex;

    m_rightFace = tempEdge.m_rightFace;
    m_leftFace  = tempEdge.m_leftFace;

    m_rightHand = tempEdge.m_rightHand;
    m_leftHand  = tempEdge.m_leftHand;
    m_rightLeg  = tempEdge.m_rightLeg;
    m_leftLeg   = tempEdge.m_leftLeg;

    m_geometry = tempEdge.m_geometry;
	m_ID       = tempEdge.m_ID;
    m_visited  = tempEdge.m_visited;

    return *this;
}


template<class VertexGeom, class EdgeGeom, class FaceGeom>
bool WEdge<VertexGeom,EdgeGeom,FaceGeom>::operator==(const WEdge<VertexGeom,EdgeGeom,FaceGeom>& tempEdge) const
{
    if( &tempEdge == this ) {
        return true;
    }

    if( m_rightHand == tempEdge.getRightHand() &&
        m_leftHand  == tempEdge.getLeftHand() &&
        m_rightLeg  == tempEdge.getRightLeg() &&
        m_leftLeg   == tempEdge.getLeftLeg() &&
        m_rightFace == tempEdge.getRightFace() &&
        m_leftFace  == tempEdge.getLeftFace() &&
        m_startVertex == tempEdge.getStartVertex() &&
        m_endVertex   == tempEdge.getEndVertex() &&
		m_ID == tempEdge.getID() )
    {
        return true;
    }
    else
    {
        return false;
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WVertex<VertexGeom,EdgeGeom,FaceGeom>* WEdge<VertexGeom,EdgeGeom,FaceGeom>::getVertexOfLeftHand() const
{
    if ( m_endVertex == m_leftHand->m_startVertex ) {
        return m_leftHand->m_endVertex;
    }
    else {
        return m_leftHand->m_startVertex;
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WVertex<VertexGeom,EdgeGeom,FaceGeom>* WEdge<VertexGeom,EdgeGeom,FaceGeom>::getVertexOfRightHand() const
{
    if ( m_endVertex == m_rightHand->m_startVertex ) {
        return m_rightHand->m_endVertex;
    }
    else {
        return m_rightHand->m_startVertex;
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::connectRightHand(WEdge<VertexGeom,EdgeGeom,FaceGeom>* rightHand)
{
    m_rightHand = rightHand;

    if ( m_endVertex == rightHand->m_startVertex ) {
        rightHand->m_rightLeg = this;
    }
    else {
        rightHand->m_leftHand = this;
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::connectLeftHand( WEdge<VertexGeom,EdgeGeom,FaceGeom>* leftHand)
{
    m_leftHand = leftHand;

    if ( m_endVertex == leftHand->m_startVertex ) {
        leftHand->m_leftLeg   = this;
    }
    else {
        leftHand->m_rightHand = this;
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::connectRightLeg( WEdge<VertexGeom,EdgeGeom,FaceGeom>* rightLeg)
{
    m_rightLeg = rightLeg;

    if ( m_startVertex == rightLeg->m_startVertex ) {
        rightLeg->m_leftLeg = this;
    }
    else {
        rightLeg->m_rightHand = this;
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::connectLeftLeg(  WEdge<VertexGeom,EdgeGeom,FaceGeom>* leftLeg)
{
    m_leftLeg = leftLeg;

    if ( m_startVertex == leftLeg->m_startVertex ) {
        leftLeg->m_rightLeg = this;
    }
    else {
        leftLeg->m_leftHand = this;
    }
}


template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::connectRightFace( WFace<VertexGeom,EdgeGeom,FaceGeom>* rightFace)
{
    m_rightFace = rightFace;

    if ( rightFace != rg_NULL ) {
        rightFace->setFirstEdge( this );
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::connectLeftFace(  WFace<VertexGeom,EdgeGeom,FaceGeom>* leftFace)
{
    m_leftFace = leftFace;

    if ( leftFace != rg_NULL ) {
        leftFace->setFirstEdge( this );
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::reverse()
{
    WVertex<VertexGeom,EdgeGeom,FaceGeom>* startVertex = m_startVertex;
    WVertex<VertexGeom,EdgeGeom,FaceGeom>* endVertex   = m_endVertex;

    WEdge<VertexGeom,EdgeGeom,FaceGeom>*   leftHand    = m_leftHand;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>*   rightHand   = m_rightHand;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>*   leftLeg     = m_leftLeg;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>*   rightLeg    = m_rightLeg;

    WFace<VertexGeom,EdgeGeom,FaceGeom>*   leftFace    = m_leftFace;
    WFace<VertexGeom,EdgeGeom,FaceGeom>*   rightFace   = m_rightFace;


    m_startVertex = endVertex;
    m_endVertex   = startVertex;

    m_leftHand    = rightLeg;
    m_rightHand   = leftLeg;
    m_leftLeg     = rightHand;
    m_rightLeg    = leftHand;

    m_leftFace    = rightFace;
    m_rightFace   = leftFace;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::getEdgesInStar(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList)
{
    m_startVertex->getEdgesInStar( edgeList );

    rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*> edgesInStarOfEndVertex;
    m_endVertex->getEdgesInStar( edgesInStarOfEndVertex );

    edgesInStarOfEndVertex.reset4Loop();
    while ( edgesInStarOfEndVertex.setNext4Loop() ) {
        WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge = edgesInStarOfEndVertex.getEntity();

        if (    currEdge != m_leftHand 
             && currEdge != m_rightHand 
             && currEdge != m_leftLeg 
             && currEdge != m_rightLeg ) {
            edgeList.add( currEdge );
        }
    }
    //rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*> edgesInStarOfStartVertex;
    //rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*> edgesInStarOfEndVertex;

    //m_startVertex->getEdgesInStar( edgesInStarOfStartVertex );
    //m_endVertex->getEdgesInStar(   edgesInStarOfEndVertex );



    //edgeList.add( m_leftHand );
    //edgeList.add( m_rightHand );
    //edgeList.add( m_leftLeg );
    //edgeList.add( m_rightLeg );

    //edgesInStarOfStartVertex.reset4Loop();
    //while ( edgesInStarOfStartVertex.setNext4Loop() ) {
    //    WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge = edgesInStarOfStartVertex.getEntity();

    //    if (    currEdge != m_leftHand 
    //         && currEdge != m_rightHand 
    //         && currEdge != m_leftLeg 
    //         && currEdge != m_rightLeg ) {
    //        edgeList.add( currEdge );
    //    }
    //}

    //edgesInStarOfEndVertex.reset4Loop();
    //while ( edgesInStarOfEndVertex.setNext4Loop() ) {
    //    WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge = edgesInStarOfEndVertex.getEntity();

    //    if (    currEdge != m_leftHand 
    //         && currEdge != m_rightHand 
    //         && currEdge != m_leftLeg 
    //         && currEdge != m_rightLeg ) {
    //        edgeList.add( currEdge );
    //    }
    //}
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WEdge<VertexGeom,EdgeGeom,FaceGeom>::getFacesInStar(rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*>& faceList)
{
    m_startVertex->getIncidentWFaceList( faceList );

    rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*> facesIncidentToEndVertex;
    m_endVertex->getIncidentWFaceList( facesIncidentToEndVertex );

    facesIncidentToEndVertex.reset4Loop();
    while ( facesIncidentToEndVertex.setNext4Loop() ) {
        WFace<VertexGeom,EdgeGeom,FaceGeom>* currFace = facesIncidentToEndVertex.getEntity();

        if ( currFace != m_leftFace && currFace != m_rightFace ) {
            faceList.add( currFace );
        }
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WFace<VertexGeom,EdgeGeom,FaceGeom>::WFace()
{
    m_firstEdge = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
	m_ID = -1;
    m_visited = rg_FALSE;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WFace<VertexGeom,EdgeGeom,FaceGeom>::WFace(const  rg_INT&		tID)
{
	m_ID = tID;
    m_firstEdge = (WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_visited = rg_FALSE;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
WFace<VertexGeom,EdgeGeom,FaceGeom>::WFace(const  rg_INT&		tID,
                                           WEdge<VertexGeom,EdgeGeom,FaceGeom>* tFirstEdge)
{
	m_ID = tID;
    m_firstEdge=tFirstEdge;
    m_visited = rg_FALSE;
}


//오류수정 : 이창희 (m_firstEdge=tGeometry; -> m_geometry=tGeometry;)
template<class VertexGeom, class EdgeGeom, class FaceGeom>
WFace<VertexGeom,EdgeGeom,FaceGeom>::WFace(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tFirstEdge, 
                                           const FaceGeom& tGeometry, const rg_INT& tID)
{
    m_firstEdge=tFirstEdge;
    m_geometry=tGeometry;
	m_ID = tID;
    m_visited = rg_FALSE;
}

//오류수정 : 이창희 (m_firstEdge=tempFace.m_firstEdge; 중복 삭제, m_geometry=tempFace.m_geometry; 삽입)
template<class VertexGeom, class EdgeGeom, class FaceGeom>
WFace<VertexGeom,EdgeGeom,FaceGeom>::WFace(const WFace<VertexGeom,EdgeGeom,FaceGeom>& tempFace)
{
    m_firstEdge=tempFace.m_firstEdge;
    m_geometry=tempFace.m_geometry;
    //m_geometry=tGeometry;
	m_ID = tempFace.m_ID;
    m_visited = tempFace.m_visited;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
WFace<VertexGeom,EdgeGeom,FaceGeom>::~WFace()
{
}


/*
template<class VertexGeom, class EdgeGeom, class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>* WFace<VertexGeom,EdgeGeom,FaceGeom>::getFirstEdge() const
{
    return m_firstEdge;
}
*/


template<class VertexGeom, class EdgeGeom, class FaceGeom>
FaceGeom WFace<VertexGeom,EdgeGeom,FaceGeom>::getGeometry() const
{
    return m_geometry;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
FaceGeom* WFace<VertexGeom,EdgeGeom,FaceGeom>::getpGeometry() 
{
    return &m_geometry;
}


template<class VertexGeom, class EdgeGeom, class FaceGeom>
rg_INT WFace<VertexGeom,EdgeGeom,FaceGeom>::getID() const
{
	return m_ID;
}



/*
template<class VertexGeom, class EdgeGeom, class FaceGeom>
rg_BOOL WFace<VertexGeom,EdgeGeom,FaceGeom>::isVisited() const
{
    return m_visited;
}
*/



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WFace<VertexGeom,EdgeGeom,FaceGeom>::isVisited(const rg_BOOL& visited)
{
    m_visited = visited;
}




template<class VertexGeom, class EdgeGeom, class FaceGeom>
rg_BOOL WFace<VertexGeom,EdgeGeom,FaceGeom>::isThere(WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex) 
{
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            if ( vertex == currEdge->getEndVertex() ) {
                return rg_TRUE;
            }
            else {
                currEdge = currEdge->getLeftHand();
            }
        }
        else {
            if ( vertex == currEdge->getStartVertex() ) {
                return rg_TRUE;
            }
            else {            
                currEdge = currEdge->getRightLeg();
            }
        }
    } while ( startEdge != currEdge );   

    return rg_FALSE;
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
rg_BOOL WFace<VertexGeom,EdgeGeom,FaceGeom>::isIncidentTo(WEdge<VertexGeom,EdgeGeom,FaceGeom>* edge) const
{
    if ( this == edge->getLeftFace() || this == edge->getRightFace() ) {
        return rg_TRUE;
    }
    else {
        return rg_FALSE;
    }
}



template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WFace<VertexGeom,EdgeGeom,FaceGeom>::setFirstEdge(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge)
{
    m_firstEdge=tempEdge;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WFace<VertexGeom,EdgeGeom,FaceGeom>::setGeometry(const FaceGeom& temp)
{
    m_geometry=temp;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WFace<VertexGeom,EdgeGeom,FaceGeom>::setID(const rg_INT& tID)
{
	m_ID = tID;
}

template<class VertexGeom, class EdgeGeom, class FaceGeom>
WFace<VertexGeom,EdgeGeom,FaceGeom>& WFace<VertexGeom,EdgeGeom,FaceGeom>::operator=(const WFace<VertexGeom,EdgeGeom,FaceGeom>& tempFace)
{
    if(this==&tempFace)
    {
        return *this;
    }

    m_firstEdge = tempFace.m_firstEdge;
    m_geometry  = tempFace.m_geometry;
	m_ID        = tempFace.m_ID;
    m_visited   = tempFace.m_visited;

    return *this;
}

template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_FLAG WFace<VertexGeom,EdgeGeom,FaceGeom>::getIncidentVertexList(rg_dList<WVertex<VertexGeom,EdgeGeom,FaceGeom>*>& vertexList)
{
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            vertexList.add( currEdge->getEndVertex() );
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            vertexList.add( currEdge->getStartVertex() );
            currEdge = currEdge->getRightLeg();
        }
        else  {
            return rg_FALSE;
        }
    } while ( startEdge != currEdge );   

    return rg_TRUE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_FLAG WFace<VertexGeom,EdgeGeom,FaceGeom>::getIncidentEdgeList(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList)
{
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do  {
        edgeList.add( currEdge );

        if ( this == currEdge->getLeftFace() )
            currEdge = currEdge->getLeftHand();
        else if ( this == currEdge->getRightFace() )  
            currEdge = currEdge->getRightLeg();
        else
            return rg_FALSE;

    } while ( startEdge != currEdge );   

    return rg_TRUE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_FLAG WFace<VertexGeom,EdgeGeom,FaceGeom>::getIncidentFaceList(rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*>& faceList)
{
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            faceList.add( currEdge->getRightFace() );
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            faceList.add( currEdge->getLeftFace() );
            currEdge = currEdge->getRightLeg();
        }
        else  {
            return rg_FALSE;
        }
    } while ( startEdge != currEdge );   

    return rg_TRUE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_FLAG WFace<VertexGeom,EdgeGeom,FaceGeom>::getBoundingVertexList(rg_dList<WVertex<VertexGeom,EdgeGeom,FaceGeom>*>& boundingVertices)
{
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do  {
        if ( this == currEdge->getLeftFace() )  {
            boundingVertices.add( currEdge->getEndVertex() );
            currEdge = currEdge->getLeftHand();
        }
        else if ( this == currEdge->getRightFace() )  {
            boundingVertices.add( currEdge->getStartVertex() );
            currEdge = currEdge->getRightLeg();
        }
        else  {
            return rg_FALSE;
        }
    } while ( startEdge != currEdge );   

    return rg_TRUE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_FLAG WFace<VertexGeom,EdgeGeom,FaceGeom>::getBoundingEdgeList(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& boundingEdges)
{
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
    WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do  {
        boundingEdges.add( currEdge );

        if ( this == currEdge->getLeftFace() )
            currEdge = currEdge->getLeftHand();
        else if ( this == currEdge->getRightFace() )  
            currEdge = currEdge->getRightLeg();
        else
            return rg_FALSE;

    } while ( startEdge != currEdge );   

    return rg_TRUE;
}




template<class VertexGeom,class EdgeGeom,class FaceGeom>
WVertex<VertexGeom,EdgeGeom,FaceGeom>::WVertex()
{
    m_firstEdge=(WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
	m_ID = -1;
    m_visited = rg_FALSE;
}


template<class VertexGeom,class EdgeGeom,class FaceGeom>
WVertex<VertexGeom,EdgeGeom,FaceGeom>::WVertex(const rg_INT& tID)
{
	m_ID        = tID;
    m_firstEdge =(WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_visited   = rg_FALSE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
WVertex<VertexGeom,EdgeGeom,FaceGeom>::WVertex(const rg_INT& tID, const VertexGeom& tGeometry)
{
    m_firstEdge =(WEdge<VertexGeom,EdgeGeom,FaceGeom>*)rg_NULL;
    m_geometry  = tGeometry;
	m_ID        = tID;
    m_visited = rg_FALSE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
WVertex<VertexGeom,EdgeGeom,FaceGeom>::WVertex(const WVertex<VertexGeom,EdgeGeom,FaceGeom>& tempVertex)
{
    m_firstEdge = tempVertex.m_firstEdge;
    m_geometry  = tempVertex.m_geometry;
	m_ID        = tempVertex.m_ID;
    m_visited   = tempVertex.m_visited;
}

template<class VertexGeom,class EdgeGeom,class FaceGeom>
WVertex<VertexGeom,EdgeGeom,FaceGeom>::~WVertex()
{
}




//get functions
/*
template<class VertexGeom,class EdgeGeom,class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>* WVertex<VertexGeom,EdgeGeom,FaceGeom>::getFirstEdge() const
{
    return m_firstEdge;
}
*/



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_INT WVertex<VertexGeom,EdgeGeom,FaceGeom>::getID() const
{
	return m_ID;
}



/*
template<class VertexGeom,class EdgeGeom,class FaceGeom>
VertexGeom WVertex<VertexGeom,EdgeGeom,FaceGeom>::getGeometry() const
{
	return m_geometry;
}
*/



template<class VertexGeom,class EdgeGeom,class FaceGeom>
VertexGeom* WVertex<VertexGeom,EdgeGeom,FaceGeom>::getpGeometry() 
{
	return &m_geometry;
}



/*
template<class VertexGeom, class EdgeGeom, class FaceGeom>
rg_BOOL WVertex<VertexGeom,EdgeGeom,FaceGeom>::isVisited() const
{
    return m_visited;
}
*/


template<class VertexGeom, class EdgeGeom, class FaceGeom>
void WVertex<VertexGeom,EdgeGeom,FaceGeom>::isVisited(const rg_BOOL& visited)
{
    m_visited = visited;
}


//set functions
template<class VertexGeom,class EdgeGeom,class FaceGeom>
void WVertex<VertexGeom,EdgeGeom,FaceGeom>::setFirstEdge(WEdge<VertexGeom,EdgeGeom,FaceGeom>* tempEdge)
{
    m_firstEdge=tempEdge;
}

template<class VertexGeom,class EdgeGeom,class FaceGeom>
void WVertex<VertexGeom,EdgeGeom,FaceGeom>::setGeometry(const VertexGeom& tGeometry)
{
	m_geometry = tGeometry;
}

template<class VertexGeom,class EdgeGeom,class FaceGeom>
void WVertex<VertexGeom,EdgeGeom,FaceGeom>::setID(const rg_INT& tID)
{
	m_ID = tID;
}

//operator overloading
template<class VertexGeom,class EdgeGeom,class FaceGeom>
WVertex<VertexGeom,EdgeGeom,FaceGeom>& WVertex<VertexGeom,EdgeGeom,FaceGeom>::operator=(const WVertex<VertexGeom,EdgeGeom,FaceGeom>& tempVertex)
{
    if(this==&tempVertex)
    {
        return *this;
    }

    m_firstEdge = tempVertex.m_firstEdge;
    m_geometry  = tempVertex.m_geometry;
	m_ID        = tempVertex.m_ID;
    m_visited   = tempVertex.m_visited;

    return *this;
}

template<class VertexGeom,class EdgeGeom,class FaceGeom>
bool WVertex<VertexGeom,EdgeGeom,FaceGeom>::operator==(const WVertex<VertexGeom,EdgeGeom,FaceGeom>& tempVertex) const
{
    if(this == &tempVertex)
    {
        return true;
    }

    if( m_firstEdge == tempVertex.m_firstEdge &&
		m_ID == tempVertex.m_ID )
    {
        return true;
    }
    else
    {
        return false;
    }

}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_FLAG WVertex<VertexGeom,EdgeGeom,FaceGeom>::getNeighborWVertexList(rg_dList<WVertex<VertexGeom,EdgeGeom,FaceGeom>*>& vertexList)
{
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do {
        if ( this == currEdge->getStartVertex() ) {
            vertexList.add( currEdge->getEndVertex() );
            currEdge = currEdge->getLeftLeg();
        }
        else if ( this == currEdge->getEndVertex() ) {
            vertexList.add( currEdge->getStartVertex() );
            currEdge = currEdge->getRightHand();
        }
        else {
            return rg_FALSE;
        }
    } while (startEdge != currEdge);

    return rg_TRUE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_FLAG WVertex<VertexGeom,EdgeGeom,FaceGeom>::getIncidentWEdgeList(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList)
{
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do {
        edgeList.add( currEdge );
        if ( this == currEdge->getStartVertex() )
            currEdge = currEdge->getLeftLeg();
        else if ( this == currEdge->getEndVertex() )
            currEdge = currEdge->getRightHand();
        else
            return rg_FALSE;
    } while (startEdge != currEdge);

    return rg_TRUE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_FLAG WVertex<VertexGeom,EdgeGeom,FaceGeom>::getIncidentWFaceList(rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*>& faceList)
{
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge = m_firstEdge;
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge  = startEdge;

    do {
        if ( this == currEdge->getStartVertex() ) {
            if ( currEdge->getLeftFace() != rg_NULL ) {
                faceList.add( currEdge->getLeftFace() );
            }
            currEdge = currEdge->getLeftLeg();
        }
        else if ( this == currEdge->getEndVertex() ) {
            if ( currEdge->getRightFace() != rg_NULL ) {
                faceList.add( currEdge->getRightFace() );
            }
            currEdge = currEdge->getRightHand();
        }
        else {
            return rg_FALSE;
        }
    } while (startEdge != currEdge);

    return rg_TRUE;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
void WVertex<VertexGeom,EdgeGeom,FaceGeom>::getEdgesInStar(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList)
{
    rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*> edgesIncidentToVertex;
    getIncidentWEdgeList(edgesIncidentToVertex);

    rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*> edgesInShellOfVertex;
    getEdgesInShell(edgesInShellOfVertex);


    edgeList.append( edgesIncidentToVertex );
    edgeList.append( edgesInShellOfVertex );
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
void WVertex<VertexGeom,EdgeGeom,FaceGeom>::getEdgesInShell(rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*>& edgeList)
{
    rg_dList<WEdge<VertexGeom,EdgeGeom,FaceGeom>*> incidentEdgeList;
    getIncidentWEdgeList(incidentEdgeList);

    rg_INT numIncidentEdges = incidentEdgeList.getSize();
    rg_dNode< WEdge<VertexGeom,EdgeGeom,FaceGeom>* >* currNode = incidentEdgeList.getFirstpNode();

    for ( rg_INT i=0; i<numIncidentEdges; i++, currNode=currNode->getNext() ) {
        WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge = currNode->getEntity();
        WEdge<VertexGeom,EdgeGeom,FaceGeom>* nextEdge = currNode->getNext()->getEntity();

        WEdge<VertexGeom,EdgeGeom,FaceGeom>* coverEdge = rg_NULL;
        if ( this == currEdge->getStartVertex() ) {
            coverEdge = currEdge->getLeftHand();
        }
        else {
            coverEdge = currEdge->getRightLeg();
        }

        if ( this == nextEdge->getStartVertex() ) {
            if ( coverEdge != nextEdge->getRightHand() ) {
                coverEdge = rg_NULL;
            }
        }
        else {
            if ( coverEdge != nextEdge->getLeftLeg() ) {
                coverEdge = rg_NULL;
            }
        }

        if ( coverEdge != rg_NULL ) {
            edgeList.add( coverEdge );
        }
    }
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
WEdge<VertexGeom,EdgeGeom,FaceGeom>* WVertex<VertexGeom,EdgeGeom,FaceGeom>::findConnectingEdge(WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex)
{
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* connectingEdge = rg_NULL;
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* startEdge      = m_firstEdge;
	WEdge<VertexGeom,EdgeGeom,FaceGeom>* currEdge       = startEdge;

    do {
        if ( this == currEdge->getStartVertex() ) {
            if ( vertex == currEdge->getEndVertex() ) {
                connectingEdge = currEdge;
                break;
            }
            currEdge = currEdge->getLeftLeg();
        }
        else { //if ( this == currEdge->getEndVertex() ) {
            if ( vertex == currEdge->getStartVertex() ) {
                connectingEdge = currEdge;
                break;
            }

            currEdge = currEdge->getRightHand();
        }
    } while (startEdge != currEdge);

    return connectingEdge;
}



template<class VertexGeom,class EdgeGeom,class FaceGeom>
rg_INT WVertex<VertexGeom,EdgeGeom,FaceGeom>::findSharingFace(WVertex<VertexGeom,EdgeGeom,FaceGeom>* vertex,
                                                              rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*>& faceList)
{
    rg_dList<WFace<VertexGeom,EdgeGeom,FaceGeom>*> facesIncidentToThisVtx;
    getIncidentWFaceList(facesIncidentToThisVtx);

    facesIncidentToThisVtx.reset4Loop();
    while ( facesIncidentToThisVtx.setNext4Loop() ) {
        WFace<VertexGeom,EdgeGeom,FaceGeom>* currFace = facesIncidentToThisVtx.getEntity();

        if ( currFace->isThere( vertex ) ) {
            faceList.add( currFace );
        }
    }

    return faceList.getSize();
}


#endif
