#include "Polygon2D.h"
#include "rg_GeoFunc.h"
#include "FileOperations.h"
#include "StringFunctions.h"
#include <fstream>
#include <sstream>
#include <iomanip>

Polygon2D::Polygon2D()
{
	m_numberOfVertices = 0;
}


Polygon2D::Polygon2D(const list<rg_Point2D>& exteriorBoundaryVertices)
{
    m_numberOfVertices = 0;
    set_exterior_boundary_vertices(exteriorBoundaryVertices);
}


Polygon2D::Polygon2D(const rg_Point2D* exteriorBoundaryVertices, const rg_INT& numVertices)
{
    set_exterior_boundary_vertices(exteriorBoundaryVertices, numVertices);
}


Polygon2D::Polygon2D(const list< list<rg_Point2D> >& boundaryVertices)
{
    m_numberOfVertices = 0;
	set_boundary_vertices(boundaryVertices);
}


Polygon2D::Polygon2D(const Polygon2D& polygon)
{
	m_boundaryVertices = polygon.m_boundaryVertices;
	m_numberOfVertices = polygon.m_numberOfVertices;
	m_boundingBox      = polygon.m_boundingBox;
}


Polygon2D::~Polygon2D()
{
}


void Polygon2D::get_boundary_vertices(list<rg_Point2D>& boundaryVertices) const
{
	for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin(); 
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
		const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin(); 
             i_vertex != verticesOnCurrentShell.end(); 
             ++i_vertex)
		{
			boundaryVertices.push_back(*i_vertex);
		}
	}
}

rg_INT Polygon2D::get_exterior_boundary_vertices(list<rg_Point2D>& exteriorBoundaryVertices) const
{
    exteriorBoundaryVertices = m_boundaryVertices.front();
    return exteriorBoundaryVertices.size();
}

rg_INT Polygon2D::get_exterior_boundary_vertices(rg_Point2D*& exteriorBoundaryVertices) const
{
    list<rg_Point2D> exteriorBoundary = m_boundaryVertices.front();
    rg_INT numVertices = exteriorBoundary.size();
    exteriorBoundaryVertices = new rg_Point2D[numVertices];

    rg_INDEX i = 0;
    for (list<rg_Point2D>::const_iterator i_vertex = exteriorBoundary.begin();
        i_vertex != exteriorBoundary.end();
        ++i_vertex, ++i)
    {
        exteriorBoundaryVertices[ i ] = (*i_vertex);
    }

    return numVertices;
}


void Polygon2D::set_exterior_boundary_vertices(const list<rg_Point2D>& exteriorBoundaryVertices)
{
    if (m_numberOfVertices > 0)
    {
        m_boundaryVertices.front().clear();
        m_boundaryVertices.front() = exteriorBoundaryVertices;
    }
    else if (m_numberOfVertices == 0)
    {
        m_boundaryVertices.push_back(exteriorBoundaryVertices);
    }

    m_numberOfVertices = 0;
    m_boundingBox.reset();
    for (list<rg_Point2D>::const_iterator i_vertex = exteriorBoundaryVertices.begin();
         i_vertex != exteriorBoundaryVertices.end(); 
         ++i_vertex)
    {
        m_boundingBox.contain(*i_vertex);
        ++m_numberOfVertices;
    }
}


void Polygon2D::set_exterior_boundary_vertices(const rg_Point2D* exteriorBoundaryVertex, const rg_INT& numVertices)
{
    if (numVertices <= 0 || exteriorBoundaryVertex == rg_NULL)
    {
        return;
    }
       
    list<rg_Point2D> copiedExteriorBoundaryVertex;
   
    for (rg_INT i = 0; i < numVertices; ++i)
    {
        copiedExteriorBoundaryVertex.push_back(exteriorBoundaryVertex[i]);
    }
        
    set_exterior_boundary_vertices(copiedExteriorBoundaryVertex);
}


void Polygon2D::set_boundary_vertices(const list< list<rg_Point2D> >& boundaryVertices)
{
    clear();

	m_boundaryVertices = boundaryVertices;
	m_numberOfVertices = 0;

    for (list< list<rg_Point2D> >::iterator i_shell = m_boundaryVertices.begin();
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
		for (list<rg_Point2D>::iterator i_vertex = i_shell->begin(); 
             i_vertex != i_shell->end(); 
             ++i_vertex)
		{
			m_boundingBox.contain(*i_vertex);
			++m_numberOfVertices;
		}
	}
}


Polygon2D& Polygon2D::operator=(const Polygon2D& polygon)
{
    if (this != &polygon)
    {
        m_boundaryVertices = polygon.m_boundaryVertices;
        m_numberOfVertices = polygon.m_numberOfVertices;
        m_boundingBox      = polygon.m_boundingBox;
    }
    return *this;
}


void Polygon2D::clear()
{
    for (list< list<rg_Point2D> >::iterator i_shell = m_boundaryVertices.begin();
        i_shell != m_boundaryVertices.end();
        ++i_shell)
    {
        (*i_shell).clear();
    }

    m_boundaryVertices.clear();

    m_numberOfVertices = 0;
    m_boundingBox.reset();
}


bool Polygon2D::read_polygon_file(const std::string& filename)
{
    list<Polygon2D> polygons;
    bool bFirstPolygonIsContainer = false;
    bool validPolygons = Polygon2D::read_polygon_file(polygons, filename, bFirstPolygonIsContainer);

    if (validPolygons)
    {
        *this = polygons.front();
        return true;
    }
    else
        return false;

    /*
    const string POLYGON( "polygon" );
    const string SHELL( "shell" );

	const int    BUFFER_SIZE = 300;
	char	     lineData[BUFFER_SIZE];
    char         seps[4] = { ' ', '\t', '\n', '\0' };
	char*	     token      = NULL;
	char*	     next_token = NULL;


	ifstream fin;
	fin.open(filename.c_str());
	if (!fin.is_open()) {
		return;
	}

    int numPolygons;
    fin.getline( lineData, BUFFER_SIZE );
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
    token = strtok_s( lineData, seps, &next_token );
#else
    token = strtok_r( lineData, seps, &next_token );
#endif	
    if ( token != NULL ) {
        numPolygons = atoi( token );
    }


	while (!fin.eof()) 
    //for(rg_INT i=0;i < numShells;++i)
    {
		fin.getline(lineData, BUFFER_SIZE);
        #if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
            token = strtok_s(lineData, seps, &next_token);
        #else
            token = strtok_r(lineData, seps, &next_token);
        #endif	

        if (token == NULL)
            continue;
        if ( POLYGON.compare( token ) == 0 ) {
            int numShells;
            fin.getline( lineData, BUFFER_SIZE );
#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
            token = strtok_s( lineData, seps, &next_token );
#else
            token = strtok_r( lineData, seps, &next_token );
#endif	
            if ( token != NULL ) {
                numShells = atoi( token );
            }

            if ( SHELL.compare( token ) == 0 ) {

            }

        }


		if (SHELL.compare(token) == 0)
		{
			//  add new shell
			list<rg_Point2D> shellBoundaryVertices;
			
			//  shell ID
            #if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
                token = strtok_s(NULL, seps, &next_token);
            #else 
                token = strtok_r(NULL, seps, &next_token);
            #endif	
			int shellID = -1;
			if (token != NULL) 
			{
                shellID = atoi(token);
			}

			//  number of points in a shell
			fin.getline(lineData, BUFFER_SIZE);
            #if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
                token = strtok_s(lineData, seps, &next_token);
            #else
                token = strtok_r(lineData, seps, &next_token);
            #endif	
			int numPointsOfCurrentShell = atoi(token);

			//  add each point on current shell boundary
			for (int i = 0; i < numPointsOfCurrentShell; ++i)
			{
				fin.getline(lineData, BUFFER_SIZE);

                #if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
                    token = strtok_s(lineData, seps, &next_token);
                #else
                    token = strtok_r(lineData, seps, &next_token);
                #endif	
				rg_REAL xCoord = atof(token);

                #if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64) 
                    token = strtok_s(NULL, seps, &next_token);
                #else
                    token = strtok_r(NULL, seps, &next_token);
                #endif	
				rg_REAL yCoord = atof(token);

				rg_Point2D point(xCoord, yCoord);
				shellBoundaryVertices.push_back(point);

				m_boundingBox.contain(point);

				++m_numberOfVertices;
			}
			m_boundaryVertices.push_back(shellBoundaryVertices);
		}
	}
    */
}


bool Polygon2D::write_polygon_file(const std::string & filename)
{
    list<Polygon2D> polygons;
    polygons.push_back(*this);
    return Polygon2D::write_polygon_file(polygons, filename, false);

    /*
    * 
	ofstream fout;
	fout.open(filename.c_str());

    rg_INT numShells = m_boundaryVertices.size();

    fout << numShells << endl;

	string SHELL("Shell"); // SHELL 0: outer shell
	int IDOfShell = 0;

	for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin(); 
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
		int numberOfVerticesOnCurrShell = i_shell->size();
		fout << SHELL << " " << IDOfShell++ << endl;
		fout << numberOfVerticesOnCurrShell << endl;

        const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);
        list<rg_Point2D>::const_iterator i_lastVertex  = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin(); 
             i_vertex != verticesOnCurrentShell.end(); 
             ++i_vertex)
		{
			const rg_Point2D& currentVertex = *i_vertex;
			fout << currentVertex.getX() << " " << currentVertex.getY() << endl;
		}
	}
    *
    */
}

bool Polygon2D::read_polygon_file(list<Polygon2D>& polygons, const std::string& filename, bool& bFirstPolygonIsContainer)
{










    // remove blank lines and lines with comments    
    string polygonFileInString;
    const char commentMark = '#';
    FileOperations::remove_both_comment_and_blank_lines_from_ascii_file(filename, commentMark, polygonFileInString);
    
    const string POLYGON("polygon");
    const char   LF('\n');
    const char   TAB('\t');
    const char   BLANK(' ');
    string delimiters = string(1, TAB) + string(1, BLANK);
    //string delimiters = string(&TAB) + string(&BLANK); does not work on Windows

    std::istringstream istrstream(polygonFileInString);
    
    // get number of polygons
    std::string line;
    std::getline(istrstream, line, LF);
    rg_INT numPolygons = std::stoi(line);
    if (numPolygons <= 0)
        return false;

    bFirstPolygonIsContainer = false;
    while (std::getline(istrstream, line, LF))
    {
        vector<string> tokenVector;
        StringFunctions::tokens_from_string(line, delimiters, tokenVector);
        if (tokenVector.front() == POLYGON)
        {
            rg_INT idOfPolygon = std::stoi(tokenVector[1]);
            Polygon2D polygon;
            rg_Point3D RGBOfVertex, RGBOfEdge;
            bool validPolygon;

            if (idOfPolygon >= 0)
                validPolygon = read_a_polygon(istrstream, line, polygon, RGBOfVertex, RGBOfEdge);
            else
            {
                validPolygon = read_a_container_polygon(istrstream, line, polygon, RGBOfVertex, RGBOfEdge);
                bFirstPolygonIsContainer = true;
            }

            if (!validPolygon)
                return false;
            else
                polygons.push_back(polygon);
        }
        else
        {
            return false;
        }
    }

    return true;
}

bool Polygon2D::read_a_polygon(std::istringstream& istrstream, std::string& line, Polygon2D& polygon, rg_Point3D& RGBOfVertex, rg_Point3D& RGBOfEdge)
{
    const string SHELL("shell");
    const char   TAB('\t');
    const char   BLANK(' ');
    string delimiters = string(1, TAB) + string(1, BLANK);
    //string delimiters = string(&TAB) + string(&BLANK); does not work on Windows

    const char LF('\n');

    list< list<rg_Point2D> > boundaryVerticesOfPolygon;
    // get number of shells
    std::getline(istrstream, line, LF);
    rg_INT numShells = std::stoi(line);
    
    if (numShells <= 0)
        return false;

    for (rg_INT i = 0; i < numShells; ++i)
    {
        std::getline(istrstream, line, LF);
        if (line.find_first_of(SHELL) != string::npos)
        {
            list<rg_Point2D> boundaryVerticesOfShell;
            // get number of vertices
            std::getline(istrstream, line, LF);
            rg_INT numVertices = std::stoi(line);

            if (numVertices <= 0)
                return false;

            // get boundary vertices of shell
            for (rg_INT j = 0; j < numVertices; ++j)
            {
                // vertex id, x-coord, y-coord, RGB of color for vertex, and RGB of color for edge
                std::getline(istrstream, line, LF);
                vector<string> tokenVector;
                rg_INT numTokens = StringFunctions::tokens_from_string(line, delimiters, tokenVector);

                switch (numTokens)
                {
                case 9:
                    {
                        RGBOfEdge.setX(std::stod(tokenVector[6]));
                        RGBOfEdge.setY(std::stod(tokenVector[7]));
                        RGBOfEdge.setZ(std::stod(tokenVector[8]));
                    }
                case 6:
                    {
                        RGBOfVertex.setX(std::stod(tokenVector[3]));
                        RGBOfVertex.setY(std::stod(tokenVector[4]));
                        RGBOfVertex.setZ(std::stod(tokenVector[5]));
                    }
                case 3:
                    {
                        rg_INT  id     = std::stoi(tokenVector[0]);
                        rg_REAL xCoord = std::stod(tokenVector[1]);
                        rg_REAL yCoord = std::stod(tokenVector[2]);
                        rg_Point2D vertex(xCoord, yCoord);
                        boundaryVerticesOfShell.push_back(vertex);
                    }
                    break;
                default:
                    break;
                }
                /*
                stringstream ss(line);
                std::string token;
                std::getline(ss, token, TAB); // id
                std::getline(ss, token, TAB); // x-coord
                rg_REAL xCoord = std::stod(token);
                std::getline(ss, token, TAB); // y-coord
                rg_REAL yCoord = std::stod(token);
                rg_Point2D vertex(xCoord, yCoord);
                boundaryVerticesOfShell.push_back(vertex);
                */
            }
            boundaryVerticesOfPolygon.push_back(boundaryVerticesOfShell);
        }
        else
        {
            return false;
        }
    }
    polygon.set_boundary_vertices(boundaryVerticesOfPolygon);

    return true;
}

bool Polygon2D::read_a_container_polygon(std::istringstream& istrstream, std::string& line, Polygon2D& polygon, rg_Point3D& RGBOfVertex, rg_Point3D& RGBOfEdge)
{
    const char   TAB('\t');
    const char   BLANK(' ');
    string delimiters = string(1, TAB) + string(1, BLANK);
    //string delimiters = string(&TAB) + string(&BLANK); does not work on Windows
    const char LF('\n');

    list< list<rg_Point2D> > boundaryVerticesOfPolygonContainer;
    // get number of vertices
    std::getline(istrstream, line, LF);
    rg_INT numVertices = std::stoi(line);

    if (numVertices <= 0)
        return false;

    // NOTE: container polygon consists of a shell whose orientaion is CW.
    list<rg_Point2D> boundaryVertices; 

    for (rg_INT i = 0; i < numVertices; ++i)
    {
        // vertex id, x-coord, y-coord, RGB of color for vertex, and RGB of color for edge
        std::getline(istrstream, line, LF);
        vector<string> tokenVector;
        rg_INT numTokens = StringFunctions::tokens_from_string(line, delimiters, tokenVector);

        switch (numTokens)
        {
        case 9:
        {
            RGBOfEdge.setX(std::stod(tokenVector[6]));
            RGBOfEdge.setY(std::stod(tokenVector[7]));
            RGBOfEdge.setZ(std::stod(tokenVector[8]));
        }
        case 6:
        {
            RGBOfVertex.setX(std::stod(tokenVector[3]));
            RGBOfVertex.setY(std::stod(tokenVector[4]));
            RGBOfVertex.setZ(std::stod(tokenVector[5]));
        }
        case 3:
        {
            rg_INT  id = std::stoi(tokenVector[0]);
            rg_REAL xCoord = std::stod(tokenVector[1]);
            rg_REAL yCoord = std::stod(tokenVector[2]);
            rg_Point2D vertex(xCoord, yCoord);
            boundaryVertices.push_back(vertex);
        }
        break;
        default:
            break;
        }
    }

    boundaryVerticesOfPolygonContainer.push_back(boundaryVertices);
    polygon.set_boundary_vertices(boundaryVerticesOfPolygonContainer);

    return true;
}

bool Polygon2D::write_polygon_file(list<Polygon2D>& polygons, const std::string& filename, const bool& bFirstPolygonIsContainer)
{
    ofstream fout;
    fout.open(filename.c_str());
    write_description_of_polygon_file_in_commented_lines(fout);

    rg_INT numPolygons = polygons.size();
    if (numPolygons <= 0)
        return false;

    fout << numPolygons << '\n';

    rg_INDEX polygonIndex = 0;
    if (bFirstPolygonIsContainer)
        polygonIndex = -1;

    rg_INDEX vertexIndex  = 0;

    for (auto& polygon : polygons)
    {
        rg_Point3D RGBOfVertex(0.0, 0.0, 0.0), RGBOfEdge(0.0, 0.0, 0.0);
        bool validPolygon;
        
        if (polygonIndex == -1)
            validPolygon = write_a_container_polygon(fout, polygon, polygonIndex, vertexIndex, RGBOfVertex, RGBOfEdge);            
        else
            validPolygon = write_a_polygon(fout, polygon, polygonIndex, vertexIndex, RGBOfVertex, RGBOfEdge);

        if (!validPolygon)
            return false;

        ++polygonIndex;
    }

    return true;
}

void Polygon2D::write_description_of_polygon_file_in_commented_lines(ofstream& fout)
{
    fout << "\# Polygon file format (PFF) Ver 0.1 by VDRC (Refer to \"polygon\_file\_format.pff\".)" << '\n';
    fout << "\# 1.  Number of polygons." << '\n';
    fout << "\# 2.  Polygon identifier (*The minus id corresponds to a container polygon.)." << '\n';
    fout << "\# 3.  Number of shells (*The container polygon has only one shell in CCW.)." << '\n';
    fout << "\# 4.  Shell identifier (*The first one is outer shell in CCW. The others in CW.)." << '\n';
    fout << "\# 5.  Number of vertices on shell boundary." << '\n';
    fout << "\# 6.  Vertex ID, X-coord, Y-coord, R-, G-, B-coord of colors for vertex and edge." << '\n';
}

bool Polygon2D::write_a_polygon(ofstream& fout, Polygon2D& polygon, rg_INDEX& polygonIndex, rg_INDEX& vertexIndex, rg_Point3D& RGBOfVertex, rg_Point3D& RGBOfEdge)
{
    const string POLYGON("polygon");
    const string SHELL("shell");
    const string TAB("\t");
    const string BLANK(" ");

    rg_INT numShells = polygon.num_shells();
    if (numShells <= 0)
        return false;

    fout << POLYGON << BLANK << polygonIndex << '\n';
    fout << numShells << '\n';
    rg_INDEX shellIndex = 0;

    list< list<rg_Point2D> > boundaryVerticesOfShells;
    polygon.get_boundary_vertices(boundaryVerticesOfShells);

    for (auto& verticesOnShell : boundaryVerticesOfShells)
    {
        rg_INT numVerticesOnShell = verticesOnShell.size();
        if (numVerticesOnShell <= 0)
            return false;

        fout << SHELL << BLANK << shellIndex++ << '\n';
        fout << numVerticesOnShell << '\n';
        fout << std::fixed;
        for (auto& vertex : verticesOnShell)
        {
            fout << vertexIndex++        << TAB                  << std::setprecision(6)        << vertex.getX()      << TAB << std::setprecision(6) << vertex.getY()      << TAB;
            fout << std::setprecision(6) << RGBOfVertex.getX()   << TAB << std::setprecision(6) << RGBOfVertex.getY() << TAB << std::setprecision(6) << RGBOfVertex.getZ() << TAB;
            fout << std::setprecision(6) << RGBOfEdge.getX()     << TAB << std::setprecision(6) << RGBOfEdge.getZ()   << TAB << std::setprecision(6) << RGBOfEdge.getZ()   << '\n';
        }
    }

    return true;
}

bool Polygon2D::write_a_container_polygon(ofstream& fout, Polygon2D& polygon, rg_INDEX& polygonIndex, rg_INDEX& vertexIndex, rg_Point3D& RGBOfVertex, rg_Point3D& RGBOfEdge)
{
    const string POLYGON("polygon");
    const string TAB("\t");
    const string BLANK(" ");

    rg_INT numShells = polygon.num_shells();
    if (numShells <= 0)
        return false;

    fout << POLYGON << BLANK << polygonIndex << '\n';

    list<rg_Point2D> boundaryVerticesOfOuterShell;
    polygon.get_exterior_boundary_vertices(boundaryVerticesOfOuterShell);

    rg_INT numVerticesOnShell = boundaryVerticesOfOuterShell.size();
    if (numVerticesOnShell <= 0)
        return false;

    fout << numVerticesOnShell << '\n';
    fout << std::fixed;
    for (auto& vertex : boundaryVerticesOfOuterShell)
    {
        fout << vertexIndex++ << TAB << std::setprecision(6) << vertex.getX() << TAB << std::setprecision(6) << vertex.getY() << TAB;
        fout << std::setprecision(6) << RGBOfVertex.getX() << TAB << std::setprecision(6) << RGBOfVertex.getY() << TAB << std::setprecision(6) << RGBOfVertex.getZ() << TAB;
        fout << std::setprecision(6) << RGBOfEdge.getX() << TAB << std::setprecision(6) << RGBOfEdge.getZ() << TAB << std::setprecision(6) << RGBOfEdge.getZ() << '\n';
    }

    return true;
}


rg_REAL Polygon2D::find_minimum_edge_length() const
{
	rg_REAL minEdgeLength = DBL_MAX;

    for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
		const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::const_iterator i_lastVertex  = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin(); 
             i_vertex != i_lastVertex; 
             ++i_vertex)
		{
			const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
            if (currentEdgeLength < minEdgeLength)
            {
                minEdgeLength = currentEdgeLength;
            }
		}

        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();
		rg_REAL currentEdgeLength       = currentVertex.distance(nextVertex);
       
        if (currentEdgeLength < minEdgeLength)
        {
            minEdgeLength = currentEdgeLength;
        }
	}

	return minEdgeLength;

}


rg_REAL Polygon2D::compute_length_of_perimeter() const
{
    rg_REAL perimeter = 0.0;
    
    for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
    {
		const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::const_iterator i_lastVertex  = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin(); 
             i_vertex != i_lastVertex; 
             ++i_vertex)
        {
            const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));
 
            rg_REAL   currentEdgeLength = currentVertex.distance(nextVertex);
            perimeter += currentEdgeLength;
        }

        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();
        rg_REAL       currentEdgeLength = currentVertex.distance(nextVertex);

        perimeter += currentEdgeLength;
    }

    return perimeter;
}


void Polygon2D::compute_mean_and_standard_deviation_of_edge_length(rg_REAL& meanOfEdgeLength, rg_REAL& standardDeviationOfEdgeLength) const
{
	meanOfEdgeLength = compute_mean_of_edge_length();
	standardDeviationOfEdgeLength = compute_standard_deviation_of_edge_length(meanOfEdgeLength);
}


rg_REAL Polygon2D::compute_mean_of_edge_length() const
{
	return compute_mean_of_edge_length_which_is_longer_than(0.0);
	/*
	rg_REAL meanOfEdgeLength = 0.0;

	list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();

	for (; i_shell != m_boundaryVertices.end(); ++i_shell)
	{
		list<rg_Point2D> verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::iterator i_vertex = verticesOnCurrentShell.begin();
		list<rg_Point2D>::iterator i_lastVertex = --verticesOnCurrentShell.end();
		for (; i_vertex != i_lastVertex; ++i_vertex)
		{
			rg_Point2D currentVertex = *i_vertex;
			rg_Point2D nextVertex = *(++i_vertex);
			--i_vertex;

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
			meanOfEdgeLength += currentEdgeLength;
		}
		rg_Point2D currentVertex = *i_vertex;
		rg_Point2D nextVertex = *verticesOnCurrentShell.begin();
		rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
		meanOfEdgeLength += currentEdgeLength;
	}

	meanOfEdgeLength = meanOfEdgeLength / m_numberOfVertices;

	return meanOfEdgeLength;
	*/
}


rg_REAL Polygon2D::compute_standard_deviation_of_edge_length(const rg_REAL& meanOfEdgeLength) const
{
	return compute_standard_deviation_of_edge_length_which_is_longer_than(0.0, meanOfEdgeLength);
	/*
	rg_REAL standardDeviationOfEdgeLength = 0.0;

	list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();

	for (; i_shell != m_boundaryVertices.end(); ++i_shell)
	{
		list<rg_Point2D> verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::iterator i_vertex = verticesOnCurrentShell.begin();
		list<rg_Point2D>::iterator i_lastVertex = --verticesOnCurrentShell.end();
		for (; i_vertex != i_lastVertex; ++i_vertex)
		{
			rg_Point2D currentVertex = *i_vertex;
			rg_Point2D nextVertex = *(++i_vertex);
			--i_vertex;

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
			rg_REAL currentDeviation = (currentEdgeLength - meanOfEdgeLength);
			standardDeviationOfEdgeLength += (currentDeviation*currentDeviation);
		}
		rg_Point2D currentVertex = *i_vertex;
		rg_Point2D nextVertex = *verticesOnCurrentShell.begin();
		rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
		rg_REAL currentDeviation = (currentEdgeLength - meanOfEdgeLength);
		standardDeviationOfEdgeLength += (currentDeviation*currentDeviation);
	}

	standardDeviationOfEdgeLength = sqrt(standardDeviationOfEdgeLength / m_numberOfVertices);

	return standardDeviationOfEdgeLength;
	*/
}


void Polygon2D::compute_mean_and_standard_deviation_of_edge_length_which_is_longer_than(const rg_REAL & miminumLengthToMeasure, rg_REAL & meanOfEdgeLength, rg_REAL & standardDeviationOfEdgeLength) const
{
	meanOfEdgeLength = compute_mean_of_edge_length_which_is_longer_than(miminumLengthToMeasure);
	standardDeviationOfEdgeLength = compute_standard_deviation_of_edge_length_which_is_longer_than(miminumLengthToMeasure, meanOfEdgeLength);
}


rg_REAL Polygon2D::compute_mean_of_edge_length_which_is_longer_than(const rg_REAL & miminumLengthToMeasure) const
{
	rg_REAL meanOfEdgeLength = 0.0;

    for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
		const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::const_iterator i_lastVertex  = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin(); 
             i_vertex != i_lastVertex; 
             ++i_vertex)
		{
            const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
            if (currentEdgeLength > miminumLengthToMeasure)
            {
                meanOfEdgeLength += currentEdgeLength;
            }
		}

        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();
		rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
        if (currentEdgeLength > miminumLengthToMeasure)
        {
            meanOfEdgeLength += currentEdgeLength;
        }
	}

	meanOfEdgeLength = meanOfEdgeLength / m_numberOfVertices;

	return meanOfEdgeLength;
}


rg_REAL Polygon2D::compute_standard_deviation_of_edge_length_which_is_longer_than(const rg_REAL & miminumLengthToMeasure, const rg_REAL & meanOfEdgeLength) const
{
	rg_REAL standardDeviationOfEdgeLength = 0.0;

    for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
		const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::const_iterator i_lastVertex  = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin(); 
             i_vertex != i_lastVertex; 
             ++i_vertex)
		{
			const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
			if (currentEdgeLength > miminumLengthToMeasure)
			{
				rg_REAL currentDeviation = (currentEdgeLength - meanOfEdgeLength);
				standardDeviationOfEdgeLength += (currentDeviation*currentDeviation);
			}
		}

        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();

		rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
		if (currentEdgeLength > miminumLengthToMeasure)
		{
			rg_REAL currentDeviation = (currentEdgeLength - meanOfEdgeLength);
			standardDeviationOfEdgeLength += (currentDeviation*currentDeviation);
		}
	}

	standardDeviationOfEdgeLength = sqrt(standardDeviationOfEdgeLength / m_numberOfVertices);

	return standardDeviationOfEdgeLength;
}


void Polygon2D::compute_mean_and_standard_deviation_of_edge_length_which_is_shorter_than(const rg_REAL & maximumLengthToMeasure, rg_REAL & meanOfEdgeLength, rg_REAL & standardDeviationOfEdgeLength) const
{
	meanOfEdgeLength = compute_mean_of_edge_length_which_is_shorter_than(maximumLengthToMeasure);
	standardDeviationOfEdgeLength = compute_standard_deviation_of_edge_length_which_is_shorter_than(maximumLengthToMeasure, meanOfEdgeLength);
}

rg_REAL Polygon2D::compute_mean_of_edge_length_which_is_shorter_than(const rg_REAL & maximumLengthToMeasure) const
{
	rg_REAL meanOfEdgeLength = 0.0;

    for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
		const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::const_iterator i_lastVertex  = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin(); 
             i_vertex != i_lastVertex; 
             ++i_vertex)
		{
			const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
            if (currentEdgeLength < maximumLengthToMeasure)
            {
                meanOfEdgeLength += currentEdgeLength;
            }
		}

        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();

        rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
        if (currentEdgeLength < maximumLengthToMeasure)
        {
            meanOfEdgeLength += currentEdgeLength;
        }
	}

	meanOfEdgeLength = meanOfEdgeLength / m_numberOfVertices;

	return meanOfEdgeLength;
}

rg_REAL Polygon2D::compute_standard_deviation_of_edge_length_which_is_shorter_than(const rg_REAL & maximumLengthToMeasure, const rg_REAL & meanOfEdgeLength) const
{
	rg_REAL standardDeviationOfEdgeLength = 0.0;

    for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
		const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::const_iterator i_lastVertex  = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin(); 
             i_vertex != i_lastVertex; 
             ++i_vertex)
		{
            const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
			if (currentEdgeLength < maximumLengthToMeasure)
			{
				rg_REAL currentDeviation      = (currentEdgeLength - meanOfEdgeLength);
				standardDeviationOfEdgeLength += (currentDeviation*currentDeviation);
			}
		}

        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();

		rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
		if (currentEdgeLength < maximumLengthToMeasure)
		{
			rg_REAL currentDeviation = (currentEdgeLength - meanOfEdgeLength);
			standardDeviationOfEdgeLength += (currentDeviation*currentDeviation);
		}
	}

	standardDeviationOfEdgeLength = sqrt(standardDeviationOfEdgeLength / m_numberOfVertices);

	return standardDeviationOfEdgeLength;
}



void Polygon2D::split_all_edges_which_are_longer_than(const rg_REAL& minEdgeLengthToSplit, const rg_REAL& edgeLengthToSurvive)
{
    for (list< list<rg_Point2D> >::iterator i_shell = m_boundaryVertices.begin();
        i_shell != m_boundaryVertices.end();
        ++i_shell)
	{
		list<rg_Point2D>& verticesOnCurrentShell      = (*i_shell);
		list<rg_Point2D>::const_iterator i_lastVertex = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin();
             i_vertex != i_lastVertex; 
             ++i_vertex)
		{
            const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
			if (currentEdgeLength > minEdgeLengthToSplit)
			{
				rg_Point2D dirVec = (nextVertex - currentVertex).getUnitVector();
				rg_INT numSplits  = (rg_INT)(currentEdgeLength / edgeLengthToSurvive) - 1;
				for (rg_INT i = 1; i <= numSplits; i++)
				{
					rg_Point2D newVertex = currentVertex + i*edgeLengthToSurvive*dirVec;
					verticesOnCurrentShell.insert(i_vertex, newVertex);
					++m_numberOfVertices;
				}
			}
		}

        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();

		rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
		if (currentEdgeLength > minEdgeLengthToSplit)
		{
			rg_Point2D dirVec = (nextVertex - currentVertex).getUnitVector();
			rg_INT numSplits = (rg_INT)(currentEdgeLength / minEdgeLengthToSplit) - 1;
			for (rg_INT i = 1; i <= numSplits; i++)
			{
				rg_Point2D newVertex = currentVertex + i*minEdgeLengthToSplit*dirVec;
				verticesOnCurrentShell.push_back(newVertex);
				++m_numberOfVertices;
			}
		}
	}
}

void Polygon2D::split_all_edges_which_are_longer_than(const rg_REAL& minEdgeLengthToSplit)
{
    for (list< list<rg_Point2D> >::iterator i_shell = m_boundaryVertices.begin();
        i_shell != m_boundaryVertices.end();
        ++i_shell)
	{
		list<rg_Point2D>& verticesOnCurrentShell      = (*i_shell);
		list<rg_Point2D>::const_iterator i_lastVertex = prev(verticesOnCurrentShell.end());

		for (list<rg_Point2D>::const_iterator i_vertex = verticesOnCurrentShell.begin();
             i_vertex != i_lastVertex; 
             ++i_vertex)
		{
            const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));

			rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
			if (currentEdgeLength > minEdgeLengthToSplit)
			{
				rg_Point2D dirVec = (nextVertex - currentVertex).getUnitVector();
				rg_INT numInsertedVertices = (rg_INT)(currentEdgeLength / minEdgeLengthToSplit);
				rg_REAL lengthOfEdge = currentEdgeLength / (rg_REAL)(numInsertedVertices + 1);
				for (rg_INT i = 1; i <= numInsertedVertices; i++)
				{
					rg_Point2D newVertex = currentVertex + i*lengthOfEdge*dirVec;
					verticesOnCurrentShell.insert(i_vertex, newVertex);
					++m_numberOfVertices;
				}
			}
		}

        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();

		rg_REAL currentEdgeLength = currentVertex.distance(nextVertex);
		if (currentEdgeLength > minEdgeLengthToSplit)
		{
			rg_Point2D dirVec = (nextVertex - currentVertex).getUnitVector();
			rg_INT numInsertedVertices = (rg_INT)(currentEdgeLength / minEdgeLengthToSplit);
			rg_REAL lengthOfEdge = currentEdgeLength / (rg_REAL)(numInsertedVertices + 1);
			for (rg_INT i = 1; i <= numInsertedVertices; i++)
			{
				rg_Point2D newVertex = currentVertex + i*lengthOfEdge*dirVec;
				verticesOnCurrentShell.push_back(newVertex);
				++m_numberOfVertices;
			}
		}
	}
}

void Polygon2D::merge_all_consecutive_edge_which_are_on_same_line()
{
    for ( list< list<rg_Point2D> >::iterator i_shell = m_boundaryVertices.begin(); i_shell != m_boundaryVertices.end(); ++i_shell ) {
        list<rg_Point2D>* shell = &( *i_shell );

        // remove vertice which have same coordinate.
        list<rg_Point2D>::iterator i_pt = shell->begin();
        while ( i_pt != shell->end() ) {
            rg_Point2D pt = *i_pt;
            
            list<rg_Point2D>::iterator i_nextPt = i_pt;
            ++i_nextPt;
            rg_Point2D nextPt;
            if ( i_nextPt == shell->end() ) {
                nextPt = shell->front();
            }
            else {
                nextPt = *i_nextPt;
            }

            if ( pt == nextPt )
            {
                i_pt = shell->erase( i_pt );
                --m_numberOfVertices;
            }
            else
                ++i_pt;
        }

        // remove vertices which are on middle of colinear lines.
        i_pt = shell->begin();
        while ( i_pt != shell->end() ) {
            rg_Point2D pt = *i_pt;

            rg_Point2D prevPt;
            rg_Point2D nextPt;
            if ( i_pt == shell->begin() ) {
                prevPt = shell->back();
            }
            else {
                list<rg_Point2D>::iterator i_prevPt = i_pt;
                --i_prevPt;
                prevPt = *i_prevPt;
            }

            list<rg_Point2D>::iterator i_nextPt = i_pt;
            ++i_nextPt;
            if ( i_nextPt == shell->end() ) {
                nextPt = shell->front();
            }
            else {
                nextPt = *i_nextPt;
            }

            rg_Point2D nextDir = ( nextPt - pt ).getUnitVector();
            rg_Point2D prevDir = ( pt - prevPt ).getUnitVector();

            if ( nextDir == prevDir ) {
                i_pt = shell->erase( i_pt );
                --m_numberOfVertices;
            }
            else {
                ++i_pt;
            }
        }
    }
}



void Polygon2D::perturb_vertices_which_has_same_coordinate()
{
    for ( list< list<rg_Point2D> >::iterator i_shell = m_boundaryVertices.begin(); i_shell != m_boundaryVertices.end(); ++i_shell ) {
        list<rg_Point2D>* shell = &( *i_shell );
        for ( list<rg_Point2D>::iterator i_pt = shell->begin(); i_pt != shell->end(); ++i_pt ) {
            rg_Point2D* pt = &( *i_pt );

            bool there_is_a_point_which_has_same_coordinate = false;

            list<rg_Point2D>::iterator i_nextPt = i_pt;
            ++i_nextPt;
            for ( ; i_nextPt != shell->end(); ++i_nextPt ) {
                rg_Point2D* nextPt = &( *i_nextPt );  

                if ( rg_EQ( pt->getX(), nextPt->getX() ) && rg_EQ( pt->getY(), nextPt->getY() ) ) {
                    there_is_a_point_which_has_same_coordinate = true;
                    break;
                }
            }

            if ( there_is_a_point_which_has_same_coordinate ) {
                rg_Point2D currPt = *pt;
                rg_Point2D prevPt;
                rg_Point2D nextPt;

                if ( i_pt == shell->begin() ) {
                    prevPt = shell->back();
                }
                else {
                    list<rg_Point2D>::iterator i_prev = i_pt;
                    --i_prev;
                    prevPt = *i_prev;
                }

                list<rg_Point2D>::iterator i_next = i_pt;
                ++i_next;
                if ( i_next == shell->end() ) {
                    nextPt = shell->front();
                }
                else {
                    nextPt = *i_next;
                }

                rg_Point2D vecToPrev = ( prevPt - currPt ).getUnitVector();
                rg_Point2D vecToNext = ( nextPt - currPt ).getUnitVector();
                rg_Point2D vecForPerturbation = ( vecToPrev + vecToNext ).getUnitVector() * 0.001;

                pt->setPoint( currPt + vecForPerturbation );

                continue;
            }

            list< list<rg_Point2D> >::iterator i_nextShell = i_shell;
            ++i_nextShell;
            for ( ; i_nextShell != m_boundaryVertices.end(); ++i_nextShell ) {
                list<rg_Point2D>* nextShell = &( *i_nextShell );
                for ( i_nextPt = nextShell->begin(); i_nextPt != nextShell->end(); ++i_nextPt ) {
                    rg_Point2D* nextPt = &( *i_nextPt );

                    if ( rg_EQ( pt->getX(), nextPt->getX() ) && rg_EQ( pt->getY(), nextPt->getY() ) ) {
                        there_is_a_point_which_has_same_coordinate = true;
                        break;
                    }
                }

                if ( there_is_a_point_which_has_same_coordinate ) {
                    break;
                }
            }

            if ( there_is_a_point_which_has_same_coordinate ) {
                rg_Point2D currPt = *pt;
                rg_Point2D prevPt;
                rg_Point2D nextPt;

                if ( i_pt == shell->begin() ) {
                    prevPt = shell->back();
                }
                else {
                    list<rg_Point2D>::iterator i_prev = i_pt;
                    --i_prev;
                    prevPt = *i_prev;
                }

                list<rg_Point2D>::iterator i_next = i_pt;
                ++i_next;
                if ( i_next == shell->end() ) {
                    nextPt = shell->front();
                }
                else {
                    nextPt = *i_next;
                }

                rg_Point2D vecToPrev = ( prevPt - currPt ).getUnitVector();
                rg_Point2D vecToNext = ( nextPt - currPt ).getUnitVector();
                rg_Point2D vecForPerturbation = ( vecToPrev + vecToNext ).getUnitVector() * 0.01;

                pt->setPoint( currPt + vecForPerturbation );
            }
        }

    }
}



void Polygon2D::create_polygon_as_rectangle(const rg_Point2D& minPt, const rg_Point2D& maxPt)
{
    clear();

    double minX = minPt.getX();
    double minY = minPt.getY();
    double maxX = maxPt.getX();
    double maxY = maxPt.getY();

    list<rg_Point2D> fourVertices;
    fourVertices.push_back(rg_Point2D(minX, maxY));
    fourVertices.push_back(rg_Point2D(minX, minY));
    fourVertices.push_back(rg_Point2D(maxX, minY));
    fourVertices.push_back(rg_Point2D(maxX, maxY));

    set_exterior_boundary_vertices(fourVertices);
}


void Polygon2D::translate(const rg_Point2D& translationAmount)
{
    m_boundingBox.reset();

    for (list< list<rg_Point2D> >::iterator i_shell = m_boundaryVertices.begin();
        i_shell != m_boundaryVertices.end();
        ++i_shell)
    {
        list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);

        for (list<rg_Point2D>::iterator i_vertex = verticesOnCurrentShell.begin();
            i_vertex != verticesOnCurrentShell.end();
            ++i_vertex)
        {
            rg_Point2D& currVertex = *i_vertex;
            rg_Point2D newVertex   = currVertex + translationAmount;
            
            currVertex.setPoint(newVertex);
            m_boundingBox.contain(currVertex);
        }
    }

}


bool Polygon2D::contains(const rg_Point2D& queryPoint)
{
    bool outerShell = true;
    for (auto& verticesOnShell : m_boundaryVertices)
    {
        if (outerShell)
        {
            bool inside_flag = rg_GeoFunc::point_in_polygon(queryPoint, verticesOnShell);
            if (!inside_flag)
                return false;

            outerShell = false;
        }
        else
        {
            bool outside_flag = !(rg_GeoFunc::point_in_polygon(queryPoint, verticesOnShell));
            if (!outside_flag)
                return false;
        }
    }

    return true;
}

bool Polygon2D::this_disk_is_outside_of_polygon_exterior_boundary(const rg_Circle2D& disk, const rg_REAL& tolerance) const
{
    rg_Point2D diskCenter = disk.getCenterPt();

    // Check whether disk center is either inside the polygon or on its boundary.
    const list<rg_Point2D>& verticesOnExteriorBoundary = m_boundaryVertices.front();
    list<rg_Point2D>::const_iterator i_lastVertex      = prev(verticesOnExteriorBoundary.end());

    for (list<rg_Point2D>::const_iterator i_vertex = verticesOnExteriorBoundary.begin();
         i_vertex != i_lastVertex;
         ++i_vertex)
    {
        rg_Line2D line1(*i_vertex, *next(i_vertex));
        if (rg_NNEG(line1.signed_distance(diskCenter), tolerance))
        {
            return false;
        }
    }

    rg_Line2D line1(*i_lastVertex, *verticesOnExteriorBoundary.begin());
    if (rg_NNEG(line1.signed_distance(diskCenter), tolerance)) // [CYSONG, 2020-03-06] rg_NPOS -> rg_NNEG
    {
        return false;
    }

    rg_REAL radius = disk.getRadius();

    // Test intersection between disk and each line (not line segment) which constitutes polygon boundary.
    for (list<rg_Point2D>::const_iterator i_vertex = verticesOnExteriorBoundary.begin();
         i_vertex != i_lastVertex;
         ++i_vertex)
    {
        rg_Line<rg_Point2D> line2(*i_vertex, *next(i_vertex));
        if (rg_LT(rg_ABS(rg_GeoFunc::distanceLineVsPointl(line2, diskCenter)), radius, tolerance))
        {
            return false;
        }
    }

    rg_Line<rg_Point2D> line2(*i_lastVertex, *verticesOnExteriorBoundary.begin());
    if (rg_LT(rg_ABS(rg_GeoFunc::distanceLineVsPointl(line2, diskCenter)), radius, tolerance))
    {
        return false;
    }

    return true;
}


bool Polygon2D::this_disk_intersects_or_is_outside_of_polygon_exterior_boundary(const rg_Circle2D& disk, const rg_REAL& tolerance) const
{
    rg_Point2D diskCenter = disk.getCenterPt();

    // Check whether disk center is either outside the polygon or on its boundary.
    const list<rg_Point2D>& verticesOnExteriorBoundary = m_boundaryVertices.front();
    list<rg_Point2D>::const_iterator i_lastVertex      = prev(verticesOnExteriorBoundary.end());

    for (list<rg_Point2D>::const_iterator i_vertex = verticesOnExteriorBoundary.begin();
         i_vertex != i_lastVertex;
         ++i_vertex)
    {
        rg_Line2D line1(*i_vertex, *next(i_vertex));
        if (rg_NPOS(line1.signed_distance(diskCenter), tolerance))
        {
            return true;
        }
    }

    rg_Line2D line1(*i_lastVertex, *verticesOnExteriorBoundary.begin());
    if (rg_NPOS(line1.signed_distance(diskCenter), tolerance))
    { 
        return true;
    }

    rg_REAL radius = disk.getRadius();

    // Test intersection between disk and each line (not line segment) which constitutes polygon boundary.
    for (list<rg_Point2D>::const_iterator i_vertex = verticesOnExteriorBoundary.begin();
        i_vertex != i_lastVertex;
        ++i_vertex)
    {
        rg_Line<rg_Point2D> line2(*i_vertex, *next(i_vertex));
        if (rg_LT(rg_ABS(rg_GeoFunc::distanceLineVsPointl(line2, diskCenter)), radius, tolerance))
        {
            return true;
        }
    }

    rg_Line<rg_Point2D> line2(*i_lastVertex, *verticesOnExteriorBoundary.begin());
    if (rg_LT(rg_ABS(rg_GeoFunc::distanceLineVsPointl(line2, diskCenter)), radius, tolerance))
    {
        return true;
    }
        
    return false;
}


bool Polygon2D::this_polygon_is_convex(const list<rg_Point2D>& boundaryVertices, list<rg_Point2D*>& reflexVertices) const
{
    list<rg_Point2D>::const_iterator i_lastVertex = prev(boundaryVertices.end());

    for (list<rg_Point2D>::const_iterator i_vertex = next(boundaryVertices.begin());
        i_vertex != i_lastVertex;
        ++i_vertex)
    {
        if (this_polygon_vertex_is_reflex(boundaryVertices, i_vertex)) 
        {
            reflexVertices.push_back(const_cast<rg_Point2D*>(&*i_vertex));
        }
    }

    if (reflexVertices.empty())
    {
        return true;
    }
    else
    {
        return false;
    }
}

rg_REAL Polygon2D::compute_area() const
{
	rg_Point2D anchorPoint(0.0, 0.0);
	list<rg_Point2D> boundaryVertices;
	get_boundary_vertices(boundaryVertices);
	
	for (list<rg_Point2D>::iterator i_vertex = boundaryVertices.begin(); 
         i_vertex != boundaryVertices.end(); 
         ++i_vertex)
	{
		anchorPoint += (*i_vertex);
	}
	anchorPoint = anchorPoint / m_numberOfVertices;

	rg_REAL totalAreaOfPolygons = 0.0;

    for (list< list<rg_Point2D> >::const_iterator i_shell = m_boundaryVertices.begin();
         i_shell != m_boundaryVertices.end(); 
         ++i_shell)
	{
        rg_REAL areaOfCurrPolygon = 0.0;
        const list<rg_Point2D>& verticesOnCurrentShell = (*i_shell);
		list<rg_Point2D>::const_iterator i_vertex      = verticesOnCurrentShell.begin();
		list<rg_Point2D>::const_iterator i_lastVertex  = prev(verticesOnCurrentShell.end());

		for (; i_vertex != i_lastVertex; ++i_vertex)
		{
            const rg_Point2D& currentVertex = *i_vertex;
            const rg_Point2D& nextVertex    = *(next(i_vertex));
            areaOfCurrPolygon += rg_GeoFunc::computeSignedAreaOfTriangle(anchorPoint, currentVertex, nextVertex);
		}
        const rg_Point2D& currentVertex = *i_lastVertex;
        const rg_Point2D& nextVertex    = *verticesOnCurrentShell.begin();
        areaOfCurrPolygon += rg_GeoFunc::computeSignedAreaOfTriangle(anchorPoint, currentVertex, nextVertex);

        totalAreaOfPolygons += areaOfCurrPolygon;
	}

	return totalAreaOfPolygons;
}

//bool Polygon2D::isReflexVertex(const rg_INT& vertexArrayIndex)
//{
//	rg_Point2D vecCurrent2Previous(0.0, 0.0);
//	rg_Point2D vecCurrent2Next(0.0, 0.0);
//	if (0 < vertexArrayIndex < m_numVertices - 1)
//	{
//		vecCurrent2Previous = m_boundaryVertex[vertexArrayIndex - 1] - m_boundaryVertex[vertexArrayIndex];
//		vecCurrent2Next = m_boundaryVertex[vertexArrayIndex + 1] - m_boundaryVertex[vertexArrayIndex];
//	}
//	else if (vertexArrayIndex == 0)
//	{
//		vecCurrent2Previous = m_boundaryVertex[m_numVertices - 1] - m_boundaryVertex[vertexArrayIndex];
//		vecCurrent2Next = m_boundaryVertex[vertexArrayIndex + 1] - m_boundaryVertex[vertexArrayIndex];
//	}
//	// vertexArrayIndex == m_numVertices - 1
//	else
//	{
//		vecCurrent2Previous = m_boundaryVertex[m_numVertices - 1] - m_boundaryVertex[vertexArrayIndex];
//		vecCurrent2Next = m_boundaryVertex[0] - m_boundaryVertex[vertexArrayIndex];
//	}
//
//	if (vecCurrent2Next * vecCurrent2Previous > 0)
//		return false;
//	else
//		return true;
//}

#ifdef PYVORONOI
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

using Point2D = rg_Point2D;

void init_Polygon2D(py::module& m) {
	py::class_<Polygon2D, std::shared_ptr<Polygon2D>>(m, "Polygon2D")
		.def(py::init([]() { return Polygon2D(); }))
		// .def(py::init([](const double& x, const double& y) { return Polygon2D(x, y); }))
        .def(py::init([](const list<Point2D>& exteriorBoundaryVertices) { return Polygon2D(exteriorBoundaryVertices); }))
		.def(py::init([](const Polygon2D& polygon) { return Polygon2D(polygon); }))
		;
}
#endif