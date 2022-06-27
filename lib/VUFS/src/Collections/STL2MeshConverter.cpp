#include "STL2MeshConverter.h"
#include "BucketForPoints.h"

#include <fstream>
#include <map>
using namespace std;


STL2MeshConverter::STL2MeshConverter()
{
}


STL2MeshConverter::~STL2MeshConverter()
{
}


bool    STL2MeshConverter::read(const string& stl_filename)
{
    bool    isReadingSuccessful = false;
    if ( is_available_file(stl_filename) ) {
        if ( is_ASCII_STL_file(stl_filename) ) {
            isReadingSuccessful = read_ASCII_STL_file(stl_filename);
        }
        else {
            isReadingSuccessful = read_binary_STL_file(stl_filename);
        }
    }

    if ( isReadingSuccessful ) {
        check_triangles_with_same_vertex();
    }

    return isReadingSuccessful;
}

    

bool    STL2MeshConverter::convert_to(FaceVertexMesh& faceVertexMeshModel) const
{
    BucketForPoints<FaceVertexMesh::Vertex*> bucket;
    bucket.initialize( m_triangles.size()/2, m_boundingBox );

    unsigned int faceID = 0;
    unsigned int vtxID  = 0;
    for ( list<Triangle>::const_iterator i_tri=m_triangles.begin(); i_tri!=m_triangles.end(); ++i_tri ) {

        if ( i_tri->are_vertices_with_same_coordinates() ) {
            continue;
        }

        FaceVertexMesh::Face* newFaceOfMesh = faceVertexMeshModel.create_face();
        newFaceOfMesh->setID( ++faceID );
        newFaceOfMesh->set_normal( i_tri->normal() );


        rg_Point3D point[3] = {i_tri->getPoint(0), i_tri->getPoint(1), i_tri->getPoint(2) }; 
        for ( int i=0; i<3; ++i ) {
            //FaceVertexMesh::Vertex* vtxOfMesh = meshModel.find( point[i] );
            //if ( vtxOfMesh == rg_NULL ) {
            //    vtxOfMesh = meshModel.create_vertex();
            //    vtxOfMesh->setID( ++vtxID );
            //    vtxOfMesh->set_coordinate( point[i] );
            //}
            FaceVertexMesh::Vertex* vtxOfMesh = rg_NULL;
            if ( !bucket.has_same_data( point[i], vtxOfMesh, rg_SYSTEM_RES ) ) {
                vtxOfMesh = faceVertexMeshModel.create_vertex();
                vtxOfMesh->setID( ++vtxID );
                vtxOfMesh->set_coordinate( point[i] );

                bucket.insert( point[i], vtxOfMesh );
            }


            newFaceOfMesh->set_bounding_vertex( i, vtxOfMesh );
            vtxOfMesh->add_incident_face( newFaceOfMesh );
        }
    }

    return true;
}

    
    
bool    STL2MeshConverter::convert_to(PolygonMeshModel& polygonMeshModel) const
{
    FaceVertexMesh fv_meshModel;
    convert_to( fv_meshModel );
    return fv_meshModel.convert_to( polygonMeshModel );    
}



bool STL2MeshConverter::make_mesh_model_from_ply_file( const string& ply_filename, PolygonMeshModel& polygonMeshModel )
{
    const int   BUFFER_SIZE = 300;
    char	    lineData[BUFFER_SIZE];
    char 	    seps[3] = { ' ', '\t', '\n' };
    char*       token = NULL;
    char*       next_token = NULL;

    string      COMMENT( "#" );
    string      VERTEX( "VTX" );
    string      FACE( "FACE" );
    string      ENDSOLID( "endsolid" );

    ifstream fin;
    fin.open( ply_filename.c_str() );


    map<int, FaceVertexMesh::Vertex*> map_ID2Vertex;
    FaceVertexMesh faceVertexMesh;

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
    while ( !fin.eof() ) {
        fin.getline( lineData, BUFFER_SIZE );   // facet or endsolid
        token = strtok_s( lineData, seps, &next_token );

        if ( token == NULL ) {
            break;
        }

        if ( COMMENT.compare( token ) == 0 ) {
            continue;
        }
        else if ( VERTEX.compare( token ) == 0 ) {
            int     ID      = atoi( strtok_s( NULL, seps, &next_token ) );
            float   x       = atof( strtok_s( NULL, seps, &next_token ) );
            float   y       = atof( strtok_s( NULL, seps, &next_token ) );
            float   z       = atof( strtok_s( NULL, seps, &next_token ) );

            FaceVertexMesh::Vertex* vtx = faceVertexMesh.create_vertex();
            vtx->setID( ID );
            vtx->set_coordinate( rg_Point3D( x, y, z ) );
            map_ID2Vertex[ID] = vtx;
        }
        else if ( FACE.compare( token ) == 0 ) {
            int     ID      = atoi( strtok_s( NULL, seps, &next_token ) );
            int     ID_vtx1 = atoi( strtok_s( NULL, seps, &next_token ) );
            int     ID_vtx2 = atoi( strtok_s( NULL, seps, &next_token ) );
            int     ID_vtx3 = atoi( strtok_s( NULL, seps, &next_token ) );

            FaceVertexMesh::Face*   face = faceVertexMesh.create_face();
            FaceVertexMesh::Vertex* vtx1 = map_ID2Vertex[ID_vtx1];
            FaceVertexMesh::Vertex* vtx2 = map_ID2Vertex[ID_vtx2];
            FaceVertexMesh::Vertex* vtx3 = map_ID2Vertex[ID_vtx3];
            face->setID( ID );
            face->set_bounding_vertex( 0, vtx1 );
            face->set_bounding_vertex( 1, vtx2 );
            face->set_bounding_vertex( 2, vtx3 );
            vtx1->add_incident_face( face );
            vtx2->add_incident_face( face );
            vtx3->add_incident_face( face );
        }
    }
    fin.close();
    faceVertexMesh.convert_to( polygonMeshModel );

    return true;
#else
    while ( !fin.eof() ) {
        fin.getline( lineData, BUFFER_SIZE );   // facet or endsolid
        token = strtok_r( lineData, seps, &next_token );

        if ( token == NULL ) {
            break;
        }

        if ( COMMENT.compare( token ) == 0 ) {
            continue;
        }
        else if ( VERTEX.compare( token ) == 0 ) {
            int     ID      = atoi( strtok_r( NULL, seps, &next_token ) );
            float   x       = atof( strtok_r( NULL, seps, &next_token ) );
            float   y       = atof( strtok_r( NULL, seps, &next_token ) );
            float   z       = atof( strtok_r( NULL, seps, &next_token ) );

            FaceVertexMesh::Vertex* vtx = faceVertexMesh.create_vertex();
            vtx->setID( ID );
            vtx->set_coordinate( rg_Point3D( x, y, z ) );
            map_ID2Vertex[ID] = vtx;
        }
        else if ( FACE.compare( token ) == 0 ) {
            int     ID      = atoi( strtok_r( NULL, seps, &next_token ) );
            int     ID_vtx1 = atoi( strtok_r( NULL, seps, &next_token ) );
            int     ID_vtx2 = atoi( strtok_r( NULL, seps, &next_token ) );
            int     ID_vtx3 = atoi( strtok_r( NULL, seps, &next_token ) );

            FaceVertexMesh::Face* face = faceVertexMesh.create_face();
            FaceVertexMesh::Vertex* vtx1 = map_ID2Vertex[ID_vtx1];
            FaceVertexMesh::Vertex* vtx2 = map_ID2Vertex[ID_vtx2];
            FaceVertexMesh::Vertex* vtx3 = map_ID2Vertex[ID_vtx3];
            face->setID( ID );
            face->set_bounding_vertex( 0, vtx1 );
            face->set_bounding_vertex( 1, vtx2 );
            face->set_bounding_vertex( 2, vtx3 );
            vtx1->add_incident_face( face );
            vtx2->add_incident_face( face );
            vtx3->add_incident_face( face );
        }
    }
    fin.close();
    faceVertexMesh.convert_to( polygonMeshModel );

    return true;
#endif
}



bool STL2MeshConverter::write(const string& stl_filename, const PolygonMeshModel& model, const bool& ASCIIOrNot)
{
    if ( ASCIIOrNot ) {
        return write_ASCII_file(stl_filename, model);
    }
    else {
        return write_binary_file(stl_filename, model);
    }
}



//bool    STL2MeshConverter::convert_to(PolygonMeshModel& meshModel) const
//{
//    return false;
//}
//


bool    STL2MeshConverter::is_available_file(const string& stl_filename)
{
    ifstream fin;
    fin.open(stl_filename.c_str(), ifstream::binary);

    bool isAvailable = false;
    if ( fin.is_open() ) {
        isAvailable = true;
    }
    fin.close();

    return isAvailable;
}



bool    STL2MeshConverter::is_ASCII_STL_file(const string& stl_filename)
{
    ifstream fin;
    fin.open(stl_filename.c_str(), ifstream::binary);

    char* buffer = new char[6];
    fin.read( buffer, 6 );

    bool   is_ASCII = false;
    string head = buffer;
    if ( head.substr(0, 5).compare( "solid" ) == 0 ) {
        is_ASCII = true;
    }    
    fin.close();

    return is_ASCII;
}



bool    STL2MeshConverter::read_ASCII_STL_file(const string& stl_filename)
{
    const int   BUFFER_SIZE = 300;
	char	    lineData[BUFFER_SIZE];
    char 	    seps[3]    = {' ', '\t', '\n' };
	char*	    token      = NULL;
	char*       next_token = NULL;
    string      ENDSOLID("endsolid");

    ifstream fin;
    fin.open(stl_filename.c_str());
    

    //solid w:\20160422_cube test_low.stl
    //    facet normal 0.827417 0.561588 0
    //    outer loop
    //        vertex 28.183 -12.2067 -18.75
    //        vertex 30.6026 -15.7715 -21.875
    //        vertex 28.183 -12.2067 -21.875
    //    endloop
    //    endfacet
    //endsolid w:\20160422_cube test_low.stl

    bool            isEndOfSTLFile = false;
    unsigned int    numTriangles = 0;
    fin.getline( lineData, BUFFER_SIZE );   // solid

 /* [CYSONG] Add preprocessor for handling WIN and LNX, respectively. */
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
    while ( !isEndOfSTLFile ) {
        fin.getline( lineData, BUFFER_SIZE );   // facet or endsolid
        token = strtok_s( lineData, seps, &next_token );

        if ( ENDSOLID.compare( token ) == 0 ) {
            isEndOfSTLFile = true;
            break;
        }

        token = strtok_s(NULL, seps, &next_token);             // normal
        float nx = atof( strtok_s(NULL, seps, &next_token) );  // x of normal
        float ny = atof( strtok_s(NULL, seps, &next_token) );  // y of normal
        float nz = atof( strtok_s(NULL, seps, &next_token) );  // z of normal
        rg_Point3D normal( nx, ny, nz );

        fin.getline( lineData, BUFFER_SIZE );   // outer loop
        
        rg_Point3D vertex[3];
        for ( int i=0; i<3; ++i ) {
            fin.getline( lineData, BUFFER_SIZE );   

            token = strtok_s(lineData, seps, &next_token); // vertex
            float x = atof(strtok_s(NULL, seps, &next_token)); // x of vertex
            float y = atof(strtok_s(NULL, seps, &next_token)); // y of vertex
            float z = atof(strtok_s(NULL, seps, &next_token)); // z of vertex

            vertex[i].setPoint( x, y, z );

            m_boundingBox.update( vertex[i] );
        }

        fin.getline( lineData, BUFFER_SIZE );   // endloop
        fin.getline( lineData, BUFFER_SIZE );   // endfacet


        m_triangles.push_back( Triangle(normal, vertex[0], vertex[1], vertex[2]) );
        ++numTriangles;
    }

    fin.close();
    return true;

#else
      while ( !isEndOfSTLFile ) {
        fin.getline( lineData, BUFFER_SIZE );   // facet or endsolid
        token = strtok_r( lineData, seps, &next_token );

        if ( ENDSOLID.compare( token ) == 0 ) {
            isEndOfSTLFile = true;
            break;
        }

        token = strtok_r(NULL, seps, &next_token);             // normal
        float nx = atof( strtok_r(NULL, seps, &next_token) );  // x of normal
        float ny = atof( strtok_r(NULL, seps, &next_token) );  // y of normal
        float nz = atof( strtok_r(NULL, seps, &next_token) );  // z of normal
        rg_Point3D normal( nx, ny, nz );

        fin.getline( lineData, BUFFER_SIZE );   // outer loop
        
        rg_Point3D vertex[3];
        for ( int i=0; i<3; ++i ) {
            fin.getline( lineData, BUFFER_SIZE );   

            token = strtok_r(lineData, seps, &next_token); // vertex
            float x = atof(strtok_r(NULL, seps, &next_token)); // x of vertex
            float y = atof(strtok_r(NULL, seps, &next_token)); // y of vertex
            float z = atof(strtok_r(NULL, seps, &next_token)); // z of vertex

            vertex[i].setPoint( x, y, z );

            m_boundingBox.update( vertex[i] );
        }

        fin.getline( lineData, BUFFER_SIZE );   // endloop
        fin.getline( lineData, BUFFER_SIZE );   // endfacet


        m_triangles.push_back( Triangle(normal, vertex[0], vertex[1], vertex[2]) );
        ++numTriangles;
    }

    fin.close();
    return true;

#endif
}



bool    STL2MeshConverter::read_binary_STL_file(const string& stl_filename)
{
    char header[80];
    char attributeByte[2];

    ifstream fin;
    fin.open(stl_filename.c_str(), ios_base::binary);

    unsigned int numTriangles;


    //bool isSystemLittleEndian = true;
    //int  one = 0x00000001;
    //if ( ((char *)&one)[0] == 0 ) {
    //    isSystemLittleEndian = false;
    //}

    //////////////////////////////////////////////////////////////////////////////////////
    //
    //  for little endian
    //
    //if ( isSystemLittleEndian ) {
        fin.read( header, 80 );
        fin.read( (char*)&numTriangles, 4 );

        for ( unsigned int countTriangle=0; countTriangle<numTriangles; ++countTriangle ) {
            float nx, ny, nz;
            fin.read( (char*)&nx, 4 );
            fin.read( (char*)&ny, 4 );
            fin.read( (char*)&nz, 4 );
            rg_Point3D normal( nx, ny, nz );

            rg_Point3D vertex[3];
            for ( int i=0; i<3; ++i ) {
                float x, y, z;
                fin.read( (char*)&x, 4 );
                fin.read( (char*)&y, 4 );
                fin.read( (char*)&z, 4 );

                vertex[i].setPoint( x, y, z );

                m_boundingBox.update( vertex[i] );
            }
            fin.read( attributeByte, 2 );

            m_triangles.push_back( Triangle(normal, vertex[0], vertex[1], vertex[2]) );
        }
    //}


    fin.close();

    return true;
}



    
void    STL2MeshConverter::check_triangles_with_same_vertex()
{
    list< Triangle >::iterator i_tri;
    for ( i_tri=m_triangles.begin(); i_tri!=m_triangles.end(); ++i_tri ) {
        Triangle* currTriangle = &(*i_tri);

        rg_Point3D point[3] = { currTriangle->getPoint(0), currTriangle->getPoint(1), currTriangle->getPoint(2) };

        if (    point[0].isEqual( point[1], rg_SYSTEM_RES ) 
             || point[1].isEqual( point[2], rg_SYSTEM_RES ) 
             || point[2].isEqual( point[0], rg_SYSTEM_RES ) ) 
        {
            m_triangles_with_same_vertex.insert( currTriangle );
        }
    }
}



AxisAlignedBox  STL2MeshConverter::get_bounding_box() const
{
    AxisAlignedBox boundingBox;

    list< Triangle >::const_iterator i_tri;
    for ( i_tri=m_triangles.begin(); i_tri!=m_triangles.end(); ++i_tri ) {
        for ( int i=0; i<3; ++i ) {
            boundingBox.update( i_tri->getPoint(i) );
        }
    }

    return boundingBox;
}




bool STL2MeshConverter::write_ASCII_file(const string& stl_filename, const PolygonMeshModel& model)
{
    ofstream fout( stl_filename.c_str() );
    fout << "solid " << stl_filename.c_str() << endl;

    list<PolygonMeshModel::Face*> all_faces;
    model.all_faces( all_faces );

    list<PolygonMeshModel::Face*>::iterator i_face;
    for ( i_face=all_faces.begin(); i_face!=all_faces.end(); ++i_face ) {
        PolygonMeshModel::Face* currFace = *i_face;

        rg_Point3D normal = currFace->normal();
        fout << "  facet normal " << normal.getX() << " " << normal.getX() << " " << normal.getX() << endl; 
        fout << "    outer loop" << endl;

        vector<PolygonMeshModel::Vertex*> vertex;
        currFace->find_bounding_vertices( vertex );
        for ( int i=0; i<3; ++i ) {
            rg_Point3D point = vertex[i]->coordinate();
            fout << "      vertex " << point.getX() << " " << point.getX() << " " << point.getX() << endl; 
        }

        fout << "    endloop" << endl;
        fout << "  endfacet" << endl;
    }
    

    fout << "endsolid " << stl_filename.c_str() << endl;
    fout.close();

    return true;
}



bool STL2MeshConverter::write_binary_file(const string& stl_filename, const PolygonMeshModel& model)
{
    return false;
}




//////////////////////////////////////////////////////////////////////////////////////
//
//  class STL2MeshConverter::Triangle
//
STL2MeshConverter::Triangle::Triangle()
{
}


    
STL2MeshConverter::Triangle::Triangle(const rg_Point3D& normal, const rg_Point3D& vertex1, const rg_Point3D& vertex2, const rg_Point3D& vertex3)
: rg_Triangle( vertex1, vertex2, vertex3 ), 
  m_normal( normal )
{
}



STL2MeshConverter::Triangle::Triangle(const STL2MeshConverter::Triangle& triangle)
: rg_Triangle( triangle ), 
m_normal( triangle.m_normal )
{
}



STL2MeshConverter::Triangle::~Triangle()
{
}



STL2MeshConverter::Triangle&  STL2MeshConverter::Triangle::operator =(const STL2MeshConverter::Triangle& triangle)
{
    if ( this != &triangle ) {
        rg_Triangle::operator = (triangle);
        m_normal = triangle.m_normal;
    }

    return *this;
}



rg_Point3D STL2MeshConverter::Triangle::normal() const
{
    return m_normal;
}


void       STL2MeshConverter::Triangle::set_normal(const rg_Point3D& normal)
{
    m_normal = normal;
}


void        STL2MeshConverter::Triangle::set(const rg_Point3D& normal, const rg_Point3D& vertex1, const rg_Point3D& vertex2, const rg_Point3D& vertex3 )
{
    setPoint(0, vertex1);
    setPoint(1, vertex2);
    setPoint(2, vertex3);
    m_normal = normal;
}



bool        STL2MeshConverter::Triangle::are_vertices_with_same_coordinates() const
{
    rg_Point3D point[3] = { rg_Triangle::getPoint(0), rg_Triangle::getPoint(1), rg_Triangle::getPoint(2) };

    if (    point[0].isEqual( point[1], rg_SYSTEM_RES ) 
         || point[1].isEqual( point[2], rg_SYSTEM_RES ) 
         || point[2].isEqual( point[0], rg_SYSTEM_RES ) ) 
    {
        return true;
    }
    else {
        return false;
    }
}



