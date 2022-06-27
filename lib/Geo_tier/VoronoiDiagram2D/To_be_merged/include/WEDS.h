#ifndef WEDS_H_
#define WEDS_H_

#include "QTEdge.h"
#include "QTFace.h"
#include "QTVertex.h"
#include "rg_Circle2D.h"

#include <list>
#include <string>
#include <fstream>
using namespace std;

namespace BULL2D {
namespace GeometryTier {

class WEDS
{
private:
    list<QTEdge>       m_edges;
    list<QTFace>       m_faces;
    list<QTVertex>     m_vertices;
    list<rg_Circle2D>  m_circles;

public:
    WEDS(void);
    ~WEDS(void);

    inline int      get_number_of_vertices() const { return m_vertices.size(); }
    inline int      get_number_of_edges()    const { return m_edges.size(); }
    inline int      get_number_of_faces()    const { return m_faces.size(); }

    void getVertices( list<QTVertex*>& vertices );
    void getEdges( list<QTEdge*>& edges );
    void getFaces( list<QTFace*>& faces );

    //inline list<QTVertex>&  getVertices() {return m_vertices;}
    //inline list<QTEdge>&    getEdges()    {return m_edges;}
    //inline list<QTFace>&    getFaces()    {return m_faces;}

    void makeWEDS( const string& fname );

private:
    bool isThisEdgeAlreadyCreated( const int& ID_first_vertex, const int& ID_second_vertex, const int& ID_mate_face, QTEdge*& createdEdge );
};

}
}

#endif