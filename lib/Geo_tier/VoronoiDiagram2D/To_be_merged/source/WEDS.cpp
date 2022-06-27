#include "WEDS.h"
using namespace BULL2D::GeometryTier;



WEDS::WEDS(void)
{
}


WEDS::~WEDS(void)
{
}



void WEDS::makeWEDS( const string& fname )
{
//     string      fileNameWithPath = (string)fname;
//     ifstream    fin;
//     fin.open(fileNameWithPath);


    ifstream fin;
    fin.open(fname);

    char        buffer[100];
    const char* seps = " \t\n[]";

    fin.getline(buffer, 100);

    int numOfVertices   = atoi(strtok(buffer, seps));
    int numOfFaces      = atoi(strtok(NULL, seps));
    fin.getline(buffer, 100);
    

    // make vertex
    // make vertex list to array

    QTVertex**  vertexArray = new QTVertex*[numOfVertices];
    int**       faceTuple  = new int*[numOfFaces];

    m_vertices.push_back( QTVertex(0) );
    QTVertex* currVertex = &(*m_vertices.rbegin());
    vertexArray[0] = currVertex;

    for( int i = 1 ; i < numOfVertices ; i++ ) {
        m_vertices.push_back( QTVertex(i) );
        currVertex = &(*m_vertices.rbegin());
        
        fin.getline(buffer, 100);
        double x    = atof(strtok(buffer, seps));
        double y    = atof(strtok(NULL, seps));
        double r    = atof(strtok(NULL, seps));

        m_circles.push_back( rg_Circle2D(x, y, r) );
        rg_Circle2D* currCircle = &(*m_circles.rbegin());

        currVertex->setCircle( currCircle );
        vertexArray[i] = currVertex;
    }
    fin.getline(buffer, 100);


    // make face tuple
    // face tuple
    for( int i_face = 0 ; i_face < numOfFaces ; i_face++ ) {
        faceTuple[i_face] = new int[6];

        fin.getline(buffer, 100);
        faceTuple[i_face][0] = atoi(strtok(buffer, seps));
        faceTuple[i_face][1] = atoi(strtok(NULL, seps));
        faceTuple[i_face][2] = atoi(strtok(NULL, seps));
    }
    fin.getline(buffer, 100);

    for( int i_face = 0 ; i_face < numOfFaces ; i_face++ ) {
        fin.getline(buffer, 100);
        faceTuple[i_face][3] = atoi(strtok(buffer, seps));
        faceTuple[i_face][4] = atoi(strtok(NULL, seps));
        faceTuple[i_face][5] = atoi(strtok(NULL, seps));
    }






    // make face
    // make vertex list to array
    int numOfEdges = 0;

    for( int i_face = 0 ; i_face < numOfFaces ; i_face++ ) {

        m_faces.push_back( QTFace(i_face) );
        QTFace* currFace = &(*m_faces.rbegin());

        int ID_first_vertex     = faceTuple[i_face][0];
        int ID_second_vertex    = faceTuple[i_face][1];
        int ID_third_vertex     = faceTuple[i_face][2];
        int ID_first_mate_face  = faceTuple[i_face][3];
        int ID_second_mate_face = faceTuple[i_face][4];
        int ID_third_mate_face  = faceTuple[i_face][5];




        // make (or get) boundary edges of current face
        // first edge
        QTEdge* firstEdge;
        if( !isThisEdgeAlreadyCreated( ID_first_vertex, ID_second_vertex, ID_third_mate_face, firstEdge ) )
        {
            m_edges.push_back( QTEdge(numOfEdges, vertexArray[ID_first_vertex], vertexArray[ID_second_vertex]) );
            firstEdge = &(*m_edges.rbegin());
            numOfEdges++;

            firstEdge->setLeftFace( currFace );
        }
        else
        {
            firstEdge->setRightFace( currFace );
        }
        vertexArray[ID_first_vertex]->setFirstEdge( firstEdge );


        // second edge
        QTEdge* secondEdge;
        if( !isThisEdgeAlreadyCreated( ID_second_vertex, ID_third_vertex, ID_first_mate_face, secondEdge ) )
        {
            m_edges.push_back( QTEdge(numOfEdges, vertexArray[ID_second_vertex], vertexArray[ID_third_vertex]) );
            secondEdge = &(*m_edges.rbegin());
            numOfEdges++;

            secondEdge->setLeftFace( currFace );
        }
        else
        {
            secondEdge->setRightFace( currFace );
        }
        vertexArray[ID_second_vertex]->setFirstEdge( secondEdge );


        // third edge
        QTEdge* thirdEdge;
        if( !isThisEdgeAlreadyCreated( ID_third_vertex, ID_first_vertex, ID_second_mate_face, thirdEdge ) )
        {
            m_edges.push_back( QTEdge(numOfEdges, vertexArray[ID_third_vertex], vertexArray[ID_first_vertex]) );
            thirdEdge = &(*m_edges.rbegin());
            numOfEdges++;

            thirdEdge->setLeftFace( currFace );
        }
        else
        {
            thirdEdge->setRightFace( currFace );
        }
        vertexArray[ID_third_vertex]->setFirstEdge( thirdEdge );



        // BOOKKEEPING
        currFace->setFirstEdge( firstEdge );


        if( currFace == firstEdge->getLeftFace() )
        {
            firstEdge->setLeftHand( secondEdge );
            firstEdge->setLeftLeg( thirdEdge );
        }
        else
        {
            firstEdge->setRightLeg( secondEdge );
            firstEdge->setRightHand( thirdEdge );
        }


        if( currFace == secondEdge->getLeftFace() )
        {
            secondEdge->setLeftHand( thirdEdge );
            secondEdge->setLeftLeg( firstEdge );
        }
        else
        {
            secondEdge->setRightLeg( thirdEdge );
            secondEdge->setRightHand( firstEdge );
        }


        if( currFace == thirdEdge->getLeftFace() )
        {
            thirdEdge->setLeftHand( firstEdge );
            thirdEdge->setLeftLeg( secondEdge );
        }
        else
        {
            thirdEdge->setRightLeg( firstEdge );
            thirdEdge->setRightHand( secondEdge );
        }
    }




    if( vertexArray != NULL )
    {
        delete[] vertexArray;
    }
    if( faceTuple != NULL ) {
        for( int i = 0 ; i < numOfFaces ; i++ ) {
            delete[] faceTuple[i];
        }
        delete[] faceTuple;
    }

}



void WEDS::getVertices( list<QTVertex*>& vertices )
{
    list<QTVertex>::iterator i_vertex = m_vertices.begin();
    for( ; i_vertex != m_vertices.end() ; i_vertex++ ) {
        QTVertex* currVertex = &(*i_vertex);
        vertices.push_back( currVertex );
    }
}



void WEDS::getEdges( list<QTEdge*>& edges )
{
    list<QTEdge>::iterator i_edge = m_edges.begin();
    for( ; i_edge != m_edges.end() ; i_edge++ ) {
        QTEdge* currEdge = &(*i_edge);
        edges.push_back( currEdge );
    }
}



void WEDS::getFaces( list<QTFace*>& faces )
{
    list<QTFace>::iterator i_face = m_faces.begin();
    for( ; i_face != m_faces.end() ; i_face++ ) {
        QTFace* currFace = &(*i_face);
        faces.push_back( currFace );
    }
}



bool WEDS::isThisEdgeAlreadyCreated( const int& ID_first_vertex, const int& ID_second_vertex, const int& ID_mate_face, QTEdge*& createdEdge )
{
    bool isExist = false;
    QTEdge* currEdge = NULL;
    createdEdge = NULL;

    list<QTEdge>::iterator i_edge = m_edges.begin();
    for( ; i_edge != m_edges.end() ; i_edge++ ) {
        currEdge = &(*i_edge);

        if( ((currEdge->getStartVertex()->getID() == ID_first_vertex && currEdge->getEndVertex()->getID() == ID_second_vertex) ||
            (currEdge->getStartVertex()->getID() == ID_second_vertex && currEdge->getEndVertex()->getID() == ID_first_vertex)) &&
            currEdge->getLeftFace()->getID() == ID_mate_face )
        {
            isExist     = true;
            createdEdge = currEdge;
            break;
        }
    }

    return isExist;
}


