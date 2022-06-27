/**
 * \mainpage Introduction
 *
 *
 *
 * 
 * \section Brief Introduction
 *
 * Voronoi diagrams have been known for many important applications from science and engineering.
 * While the properties, algorithms, and software for the ordinary Voronoi diagrams of points have been
 * well-known, their counterparts of disks in 2D and spheres in 3D have not been sufficiently studied
 * and developed. In this paper, we present <b>V</b> which is the geometric library for Voronoi diagrams
 * and their dual structures of 2D disks and 3D spheres. Application programmers can easily develop
 * programs to construct the Voronoi diagrams and their dual structures by simply calling API-functions
 * of this library. Furthermore, the architecture of <b>V</b> is carefully designed so that application
 * programs can be completely independent of future modifications and improvements. <b>V</b> is 
 * implemented in the standard C++ language and will be freely available from Voronoi Diagram Research
 * Center at Hanyang University (http://voronoi.hanyang.ac.kr).
 *
 *
 *
 * \section Examples
 *
 * The following application program shows an example of the use of <b>V</b> library.
 * B0 (lines 0~8) includes the header file for <b>V</b> library
 * and makes all names from namespace V visible in the scope of the program with <tt>using namespace V;</tt>.
 * B1 (lines 16~22) loads a PDB model. 
 * The \c FileIO::load_PDB(...) command loads a PDB file, \e filename, transforms atoms in the PDB file into hard spheres, and stores the list of \c Sphere3d objects.
 * B2 (lines 25~31) computes the Voronoi diagram of atoms 
 * and counts the number of unbounded and bounded Voronoi cells.
 * B4 (lines 34~42) transforms the Voronoi diagram into its dual structure, the quasi-triangulation
 * and counts the number of vertices, edges, and faces on the squeezed-hull.
 * B5 (lines 45~54) extracts the beta-complex with the specific beta-value from the quasi-triangulation
 * and counts the number of vertices, edges, and faces on the boundary of the beta-complex.
 *
 * File: \c v_prog.cpp
\code{.cpp}
00  // B0 : include the header file for V library.
01  #include "FileIO.h"
02  #include "Sphere3d.h"
03  #include "AtomSetVoronoiDiagram.h"
04  #include "QuasiTriangulation.h"
05  #include "BetaComplex.h"
06  using namespace V;
07  
08  #include <iostream>
09  using namespace std;
10  
11  
12  int main(int argc, char* argv[])
13  {
14      string filename     = argv[1];
15      cout << "PDB file : " << filename << endl << endl;
16  
17      // B1 : load a PDB model. 
18      list<Sphere3d> atomsInMolecule;
19      FileIO::load_PDB(filename, atomsInMolecule);
20      if ( atomsInMolecule.size() == 0 ) {
21          return 0;
22      }
23      cout << "# atoms  : " << atomsInMolecule.size() << endl << endl;
24  
25  
26      // B2 : construct the Voronoi diagram, VD, of atoms.
27      AtomSetVoronoiDiagram VD;
28      VD.construct( atomsInMolecule );
29      cout << "In the Voronoi diagram," << endl;
30      cout << "there are ";
31      cout << VD.number_of_unbounded_VCells() << " unbounded cells and ";
32      cout << VD.number_of_bounded_VCells()   << " bounded cells." << endl << endl;
33  
34  
35      // B3 : transform VD into the quasi-triangulation, QT. 
36      QuasiTriangulation QT;
37      QT.transform( VD );
38      cout << "In the quasi-triangulation," << endl;
39      cout << "there are ";
40      cout << QT.number_of_QTVertices_on_squeezed_hull() << " vertices, ";
41      cout << QT.number_of_QTEdges_on_squeezed_hull() << " edges, and ";
42      cout << QT.number_of_QTFaces_on_squeezed_hull() << " faces ";
43      cout << "on the squeezed hull of QT." << endl << endl;
44  
45  
46      // B4 : extract the beta-complex from QT. 
47      double beta = 50;
48      BetaComplex BC;
49      BC.extract( QT, beta );
50      cout << "In the beta-complex (beta=" << beta << ")," << endl;
51      cout << "there are ";
52      cout << BC.number_of_BCVertices_on_boundary() << " vertices, ";
53      cout << BC.number_of_BCEdges_on_boundary()    << " edges, and ";
54      cout << BC.number_of_BCFaces_on_boundary()    << " faces ";  
55      cout << "on the boundary of the beta-complex." << endl;
56  
57      return 0;
58  }
\endcode
 *
 *  Output:
 *
\code{.cpp}
PDB file : ..\..\..\data\1C26.pdb

# atoms  : 268

In the Voronoi diagram,
there are 36 unbounded cells and 232 bounded cells.

In the quasi-triangulation,
there are 36 vertices, 102 edges, and 68 faces on the squeezed hull of QT.

In the beta-complex (beta=50),
there are 49 vertices, 141 edges, and 94 faces on the boundary of the beta-complex.
\endcode
 *
 *
 *
 *
 *
 * \subsection Ex1 Example: Input of generators in code
 *
 * The following program creates a Voronoi diagram of 15 spherical atoms.
 * The input atoms are directly inserted in the code.
 * Finally, the number of cells, faces, edges, and vertices of the Voronoi diagram are written to \c cout.
 *
 * File: \c v_prog_input1.cpp
\code{.cpp}
00  #include "AtomSetVoronoiDiagram.h"
01  using namespace V;
02  
03  #include <iostream>
04  using namespace std;
05
06  
07  int main(int argc, char* argv[])
08  {
09      list<Sphere3d> atoms;
10      atoms.push_back( Sphere3d( 0.0,  0.0,  0.0, 1.5) );
11      atoms.push_back( Sphere3d( 3.0,  0.0,  0.0, 1.0) );
12      atoms.push_back( Sphere3d(-3.0,  0.0,  0.0, 1.0) );
13      atoms.push_back( Sphere3d( 0.0,  3.0,  0.0, 1.0) );
14      atoms.push_back( Sphere3d( 0.0, -3.0,  0.0, 1.0) );
15      atoms.push_back( Sphere3d( 0.0,  0.0,  3.0, 1.0) );
16      atoms.push_back( Sphere3d( 0.0,  0.0, -3.0, 1.0) );
17      atoms.push_back( Sphere3d( 1.5,  1.5,  1.5, 0.2) );
18      atoms.push_back( Sphere3d( 1.5,  1.5, -1.5, 0.2) );
19      atoms.push_back( Sphere3d( 1.5, -1.5,  1.5, 0.2) );
20      atoms.push_back( Sphere3d( 1.5, -1.5, -1.5, 0.2) );
21      atoms.push_back( Sphere3d(-1.5,  1.5,  1.5, 0.2) );
22      atoms.push_back( Sphere3d(-1.5,  1.5, -1.5, 0.2) );
23      atoms.push_back( Sphere3d(-1.5, -1.5,  1.5, 0.2) );
24      atoms.push_back( Sphere3d(-1.5, -1.5, -1.5, 0.2) );
25  
26      AtomSetVoronoiDiagram VD;
27      VD.construct( atoms );
28  
29      cout << "The Voronoi diagram of " 
30           << VD.number_of_generators() << " atoms has " << endl
31           << VD.number_of_VCells()     << " cells, "
32           << VD.number_of_VFaces()     << " faces, "
33           << VD.number_of_VEdges()     << " edges, and "
34           << VD.number_of_VVertices()  << " vertices."  << endl;
35  }
\endcode
 * 
 *  Output:
 *
\code{.cpp}
The Voronoi diagram of 15 atoms has 
15 cells, 50 faces, 60 edges, and 48 vertices.
\endcode
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * \subsection Ex2 Example: Input of generators from a file
 *
 * The following program creates a Voronoi diagram of 15 spherical atoms.
 * The input atoms are read from a file.
 * Finally, the number of cells, faces, edges, and vertices of the Voronoi diagram are written to \c cout.
 *
 * File: \c v_prog_input2.cpp
\code{.cpp}
00  #include "AtomSetVoronoiDiagram.h"
01  #include "FileIO.h"
02  using namespace V;
03  
04  #include <iostream>
05  using namespace std;
06
07  
08  int main(int argc, char* argv[])
09  {
10      list<Sphere3d> atoms;
11      FileIO::load_atoms("15_atoms.txt", atoms);
12  
13      AtomSetVoronoiDiagram VD;
14      VD.construct( atoms );
15  
16      cout << "The Voronoi diagram of " 
17           << VD.number_of_generators() << " atoms has " << endl
18           << VD.number_of_VCells()     << " cells, "
19           << VD.number_of_VFaces()     << " faces, "
20           << VD.number_of_VEdges()     << " edges, and "
21           << VD.number_of_VVertices()  << " vertices."  << endl;
22  }
\endcode
 * 
 *  Output:
 *
\code{.cpp}
The Voronoi diagram of 15 atoms has 
15 cells, 50 faces, 60 edges, and 48 vertices.
\endcode
 * 
 * The input file which read by \c FileIO::load_atoms() has a very naked format of input spherical atoms for general applications.
 * The first command line in the file contains an integer denoting the number of input spheres in the file. 
 * Then, each following line contains the definition of each input sphere with 5 fields: \c f1 through \c f5. 
 * \c f1 denotes the input sphere id; 
 * \c f2, \c f3, and \c f4 denote the x-, y-, and z-coordinate of the sphere center; 
 * \c f5 denotes the radius of the sphere.
 * The following is the example of input file.
 *
 * File: \c 15_atoms.txt
\code{.cpp}
15
1	 0.0  0.0  0.0 1.5
2	 3.0  0.0  0.0 1.0
3	-3.0  0.0  0.0 1.0
4	 0.0  3.0  0.0 1.0
5	 0.0 -3.0  0.0 1.0
6	 0.0  0.0  3.0 1.0
7	 0.0  0.0 -3.0 1.0
8	 1.5  1.5  1.5 0.2
9	 1.5  1.5 -1.5 0.2
10	 1.5 -1.5  1.5 0.2
11	 1.5 -1.5 -1.5 0.2
12	-1.5  1.5  1.5 0.2
13	-1.5  1.5 -1.5 0.2
14	-1.5 -1.5  1.5 0.2
15	-1.5 -1.5 -1.5 0.2
\endcode
 *
 * 
 *
 *
 *
 *
 * \subsection Ex3 Example: Queries in Voronoi diagram
 *
 * The following program creates a Voronoi diagram of 15 spherical atoms
 * and answers various queries about the Voronoi diagram.
 * The program gets a bounded Voronoi cell and finds the bounding faces of the cell, 
 * takes one of bounding faces, and finds the bounding edges of the face.
 * 
 * 
 * 
 * File: \c v_prog_query.cpp
\code{.cpp}
00  #include "AtomSetVoronoiDiagram.h"
01  #include "FileIO.h"
02  using namespace V;
03  
04  #include <iostream>
05  using namespace std;
06
07  
08  int main(int argc, char* argv[])
09  {
10      list<Sphere3d> atoms;
11      FileIO::load_atoms("15_atoms.txt", atoms);
12  
13      AtomSetVoronoiDiagram VD;
14      VD.construct( atoms );
15  
16      list<VCell> boundedVCells;
17      VD.get_bounded_VCells( boundedVCells );
18  
19      VCell cell = *(boundedVCells.begin());
20      cout << "This bounded Voronoi cell has " << endl
21           << VD.number_of_bounding_VFaces(     cell ) << " bounding faces, "
22           << VD.number_of_bounding_VEdges(     cell ) << " bounding edges, and "
23           << VD.number_of_bounding_VVertices(  cell ) << " bounding vertices, "
24           << endl;
25  
26      cout << endl;
27  
28      list<VFace> boundingFaces;
29      VD.find_bounding_VFaces( cell, boundingFaces );
30  
31      VFace face = *(boundingFaces.begin());
32      cout << "This Voronoi face has " << endl
33           << VD.number_of_bounding_VEdges(     face ) << " bounding edges and "
34           << VD.number_of_bounding_VVertices(  face ) << " bounding vertices, "
35           << endl;
36  
37      cout << endl;
38  
39      list<VEdge> boundingEdges;
40      VD.find_bounding_VEdges( face, boundingEdges );
41  
42      VEdge edge = *(boundingEdges.begin());
43  
44      list<VFace> incidentFaces;
45      VD.find_incident_VFaces( edge, incidentFaces );
46  
47      cout << incidentFaces.size() << " faces are incident to the edge " << edge.ID() << "." << endl;
48      int index = 1;
49      list<VFace>::iterator i_face;
50      for ( i_face=incidentFaces.begin(); i_face!=incidentFaces.end(); ++i_face, ++index ) {
51          cout << "face " << index << ": f_" << i_face->ID() << endl;
52      }
53  }
\endcode
 * 
 *  Output:
 *
\code{.cpp}
This bounded Voronoi cell has 
14 bounding faces, 36 bounding edges, and 24 bounding vertices, 

This Voronoi face has 
8 bounding edges and 8 bounding vertices, 

3 faces are incident to the edge 1.
face 1: f_3
face 2: f_1
face 3: f_4
\endcode
 *
 *
 * 
 *
 *
 *
 *
 * 
 *
 *
 *
 *
 * \subsection Ex4 Example: Access to generators in Voronoi diagram
 *
 * After creating a Voronoi diagram, the following program gets a unbounded Voronoi cell,
 * accesses its generator, and writes the center and radius of the generator to \c cout.
 * In the next stage, the program gets a bounded Voronoi edge,
 * accesses its three generators, and writes the center and radius of each generator to \c cout.
 * In <b>V</b> library, user can access generators to define a vertex, an edge, a face, and a cell of the Voronoi diagram, respectively.
 * 
 * File: \c v_prog_gen.cpp
\code{.cpp}
00  #include "AtomSetVoronoiDiagram.h"
01  #include "FileIO.h"
02  using namespace V;
03  
04  #include <iostream>
05  using namespace std;
06
07  int main(int argc, char* argv[])
08  {
09      list<Sphere3d> atoms;
10      FileIO::load_atoms("15_atoms.txt", atoms);
11  
12      AtomSetVoronoiDiagram VD;
13      VD.construct( atoms );
14  
15      list<VCell> unboundedVCells;
16      VD.get_unbounded_VCells( unboundedVCells );
17  
18      VCell cell = *(unboundedVCells.begin());
19  
20      AtomGenerator generatorOfCell;
21      VD.generator_to_define( cell, generatorOfCell );
22  
23      cout << "Generators of the Voronoi cell " << cell.ID() << endl;
24      cout << "atom :  "
25           << " center = (" << generatorOfCell.atom().center().x() << ", "
26                            << generatorOfCell.atom().center().y() << ", "
27                            << generatorOfCell.atom().center().z() << "),\t"
28           << " radius = "  << generatorOfCell.atom().radius()     << endl;
29  
30      cout << endl;
31  
32  
33      list<VEdge> boundedEdges;
34      VD.get_bounded_VEdges( boundedEdges );
35  
36      VEdge edge = *(boundedEdges.begin());
37  
38      vector<AtomGenerator> generatorsOfEdge;
39      VD.generators_to_define( edge, generatorsOfEdge );
40  
41      cout << "Generators of the Voronoi edge " << edge.ID() << endl;
42      for ( int i=0; i<generatorsOfEdge.size(); ++i ) {
43          cout << "atom " << i << ":  "
44               << " center = (" << generatorsOfEdge[i].atom().center().x() << ", "
45                                << generatorsOfEdge[i].atom().center().y() << ", "
46                                << generatorsOfEdge[i].atom().center().z() << "),\t"
47               << " radius = "  << generatorsOfEdge[i].atom().radius()     << endl;
48      }
49  }
\endcode
 *
 *  Output:
 *
\code{.cpp}
Generators of the Voronoi cell 1
atom :   center = (3, 0, 0),	 radius = 1

Generators of the Voronoi edge 1
atom 0:   center = (0, -3, 0),	 radius = 1
atom 1:   center = (0, 0, -3),	 radius = 1
atom 2:   center = (0, 0, 0),	 radius = 1.5
\endcode
 *
 * 
 *
 *
 * \section downloanInstall Download and Installation
 *
 * This document describes how to download and install the <b>V</b> library for Windows or Linux systems.
 * Before you install the <b>V</b> library, you need to identify the correct distribution file for your platform.
 * After you have downloaded the <b>V</b> library, you have to unpack it at your installation directory.
 * The root directory, V, contains the following five subdirectories: 
 *  - \c doc : documents for the <b>V</b> library.
 *  - \c examples : sample data and C++ example programs using APIs of the <b>V</b> library.
 *  - \c include : header files of the <b>V</b> library.
 *  - \c lib : library file.
 *
 * The following figure shows the directory structure of the <b>V</b> library.
 *
 * \image html structure_of_V_installation_directory.jpg "The structure of V installation directories." width=600px
 * \image latex structure_of_V_installation_directory.pdf "The structure of V installation directories." width=10cm 
 *
 * 
 *
 *
 * \section how2useBULL How to Use V
 *
 * \subsection useBULLVS2010 Using V with Microsoft Visual C++ at Windows
 *
 * This section describes how to use the <b>V</b> library with Microsoft Visual Studio 2010 to
 *  - build and run the C++ examples.
 *  - create a C++ project and link the target with the <b>V</b> library.
 *
 * Throughout this section, <b>V</b> installation directory is refered to as <tt><VDIR></tt>.
 * 
 *
 * \subsubsection bulllibraries Libraries
 * 
 * The <b>V</b> library is precompiled and distributed in four static formats for Visual Studio 2010: 
 * multi-threaded, multi-threaded DLL, multi-threadeded debug and multi-threaded DLL debug, 
 * so that a static executable can be linked with \c libcmt.lib, \c msvcrt.lib, \c libcmtd.lib, or \c msvcrtd.lib. 
 * These four formats use the standard templete library (STL) and are compiled using the namespace \c std.
 *
 * The libraries for Visual Studio 2010 can be found in the following directories:
 *  - multi-threaded (MT)
 *      - <tt><VDIR>\\lib\\vs2017\\stat_mt\\V.lib</tt>
 *  - multi-threaded DLL (MD)
 *      - <tt><VDIR>\\lib\\vs2017\\stat_md\\V.lib</tt>
 *  - multi-threaded debug (MTD)
 *      - <tt><VDIR>\\lib\\vs2017\\stat_mtd\\V.lib</tt>
 *  - multi-threaded DLL debug (MDD)
 *      - <tt><VDIR>\\lib\\vs2017\\stat_mdd\\V.lib</tt>
 *
 *
 * \subsubsection buildExample Building and running V example
 *
 * The <b>V</b> example can be built and run for each type of static formats. 
 * The instructions below use the multi-threaded (MT) format 
 * and show the processes to build and run the example program \c v_prog.cpp for the Visual Studio 2017 environment.
 * The necessary file for the example is <tt><VDIR>\\examples\\vs2017\\v_examples\\v_examples.sln</tt>.
 * Note that the order of these instructions is important.
 *
 *  -# Start Microsoft Visual Studio 2017.
 *
 *  -# From the \b File menu, choose <b>Open Project/Solution</b>. The <b>Open Project</b> dialog box appears.
 *      - Select the directory <tt><VDIR>\\examples\\vs2017\\v_examples</tt>. 
 *      - Select the \c v_examples.sln file and click \b Open. 
 *
 *  -# From the \b Build menu, choose \b Build \b Solution. 
 *
 *      Wait until the completion of the building process. 
 *
 *  -# To run the example, press \b Ctrl+F5. The result is then displayed.
 *   
 *
 *   
 *
 *   
 *
 * \subsubsection buildYourProject Building your own project which uses V
 *
 * The following information applies to multi-threaded library for the Visual Studio 2010. 
 * Let us assume that you want to build an executable named \c myapp.exe and have: 
 *  - a source file named \c myapp.cpp which uses APIs of the <b>V</b> library, and
 *  - a directory where this file is located and which, for the sake of simplicity, we'll refer to as <tt><MYAPPDIR></tt>. 
 *
 * One recommended way to build the program is to create a project named \c myapp.sln as described below. 
 * Note that the order of instructions is important.
 *
 *  -# Start Microsoft Visual Studio 2017.
 *
 *  -# Create a solution file \c myapp.sln from the \b File menu, select <b>New > Project...</b>.
 *
 *  -# When the \b New \b Project dialog box appears, do the following steps.
 *      - In the \b Project \b types pane, select <b>Visual C++</b> and \b Win32. 
 *      - In the \b Templates pane, select the <b>Win32 Console Application</b> icon. 
 *      - Fill in the project name (\c myapp). 
 *      - If necessary, correct the location of the project (to <tt><MYAPPDIR></tt>). 
 *      - Click \b OK. 
 *
 *  -# When the <b>Win32 Application Wizard</b> appears, do the following steps.
 *      - Click on <b>Application Settings</b>. 
 *      - Select <b>Console Application</b> as <b>Application type</b>.  
 *      - Make sure that <b>Empty project</b> is checked in <b>Additional Options</b>.
 *      - Click \b Finish. 
 *
 *      Step 4 creates the solution file \c myapp.sln with a single project \c myapp. 
 *      You can see the contents of the solution by selecting <b>Solution Explorer</b> in the \b View menu.
 *
 *  -# Add your own source file to the project. From the \b Project menu, choose <b>Add Existing Item...</b>.  
 *      - Move to the directory <tt><MYAPPDIR></tt> and select \c myapp.cpp. 
 *      - Click \b Open. 
 *
 *  -# Set some options to let the project know the locations of header files and library file of <b>V</b>.
 *
 *      From the \b Project menu, choose \b Properties. 
 *      The \b myapp <b>Property Pages</b> dialog box appears.
 *
 *      In the \b Configuration drop-down list, select \b Release.
 *
 *      Select \b C/C++ in the <b>Configuration Properties</b> tree. 
 *      - Select \b General: 
 *          - In the <b>Additional Include Directories</b> field, add the directories below: 
 *              - <tt><VDIR>\\include</tt> 
 *      - Select <b>Code Generation</b>:  
 *          - Set <b>Runtime Library</b> to <b>Multi-threaded (/MT)</b>.
 *
 *      Select \b Linker in the <b>Configuration Properties</b> tree. 
 *      - Select \b General and then select <b>Additional Library Directoriess</b>. Add the directories: 
 *          - <tt><VDIR>\\lib\\vs2017\\stat_mt</tt>
 *      - Select \b Input and then select <b>Additional Dependencies</b>. Add the files: 
 *          - \c V.lib
 *
 *      Click \b OK to close the <b>myapp Property Pages</b> dialog box. 
 *
 *  -# Set the default project configuration.
 *      From the \b Build menu, select <b>Configuration Manager...</b>. 
 *      - Select \b Release in the <b>Active Solution Configuration</b> drop-down list.
 *      - Click \b Close.
 *
 *  -# Finally, to build the project, from the \b Build menu, select <b>Build Solution</b>. 
 *
 * After completion of the compiling and linking process, the executable \c myapp.exe is created at <tt><MYAPPDIR>\\myapp\\Release\\</tt>. 
 *
 *
 *
 * \subsection useBULLgcc Using V with GCC at Linux
 *
 * This section describes how to use the <b>V</b> library with GCC at Linux environment to:
 *  - build and run the C++ examples delivered with the <b>V</b> library.
 *  - create your own C++ program with the <b>V</b> library.
 *
 * Throughout this section, <b>V</b> installation directory is refered to as <tt><VDIR></tt>.
 * 
 *
 * \subsubsection libraries Libraries
 * 
 * The <b>V</b> library is provided in a static format for GCC. 
 * The libraries are precompiled by gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1 
 * and can be found in the following directories:
 *  - <tt><VDIR>\\lib\\libV.a</tt>
 *
 *
 * \subsubsection buildExampleLinux Building and running V example
 *
 * <b>V</b> examples can be built and run with the <b>V</b> library in Linux. 
 * The instructions below show the process to build and run the example code \c v_prog.cpp in Linux environment.
 * The two necessary files for the example can be found in <tt><VDIR>\\examples\\gcc\\v_prog</tt>.
 *  - \c v_prog.cpp : the example source code.
 *  - \c Makefile   : the make file to create an executable.
 *
 * You can create the executable by just running \c Makefile.
 *
 *  -# Run \c make. 
 *
 *      <tt>\$ make</tt>
 *
 *      Then, the executable \c v_prog for \c v_prog.cpp is created.
 *
 *  -# Run the executable.
 *
 *      <tt>\$ ./v_prog ../../data/1AOJ.pdb</tt>
 *
 *      The result is then displayed.
 *      The data \c 1AOJ.pdb can be found in <tt><VDIR>\\examples\\data</tt>. 
 *
 *
 *
 *
 * \subsubsection buildYourLinuxProgram Building your own program using V
 *
 * Let's assume that you want to build the executable named \c myapp and have: 
 *  - a \c C++ source file named \c myapp.cpp which uses the API of the <b>V</b> library.
 *  - a directory where this file is located and which, for the sake of simplicity, we'll refer to as <tt><MYAPPDIR></tt>. 
 *
 * One way to achieve that goal is to generate a \c Makefile at <tt><MYAPPDIR></tt>. 
 * The following \c Makefile is for \c myapp.cpp at Linux environment.
 *
 * File: \c Makefile
 *
\code{.cpp}
00  #---------------------------------------------------------------------#
01  # Makefile for myapp.cpp
02  #---------------------------------------------------------------------#
03
04  #---------------------------------------------------------------------#
05  # 1 Prepare the environment for making the executable
06  #
07  #    1.1 Set directory information
08  #---------------------------------------------------------------------#
09  
10  V_INCDIR = <VDIR>/include
11
12  
13  V_LIBDIR = <VDIR>/lib
14
15  
16  #---------------------------------------------------------------------#
17  #    1.2 Set compiler flags
18  #---------------------------------------------------------------------#
19  
20  CXX         = g++  
21  CXXFLAGS32  = -Wall -O2 -m32
22  INCLUDE_DIR = -I$(V_INCDIR)  
23  
24  #---------------------------------------------------------------------#
25  #    1.3 Set linker flags
26  #---------------------------------------------------------------------#
27  
28  LDFLAGS32 = \
29  	-Xlinker --start-group \
30  		-L$(V_LIBDIR) -lV \
31  	-Xlinker --end-group
32
33  	
34  #---------------------------------------------------------------------#
35  #    1.4. Define source files 
36  #---------------------------------------------------------------------#
37  
38  APP_SRC  = myapp.cpp
39  TARGET32 = myapp
40  
41  
42  #---------------------------------------------------------------------#
43  # 2. Create  
44  #---------------------------------------------------------------------#
45  
46  all: $(TARGET32) 
47  
48  .PHONY : $(TARGET32)
49  $(TARGET32): $(TARGET32).o
50  	${CXX} $(CXXFLAGS32) $(TARGET32).o $(LDFLAGS32) -o $(TARGET32)
51  
52  $(TARGET32).o: $(APP_SRC)
53  	${CXX} -c $(CXXFLAGS32) $(APP_SRC) $(INCLUDE_DIR) -o $(TARGET32).o
54  
55  
56  #---------------------------------------------------------------------#
57  # 3. Clean 
58  #---------------------------------------------------------------------#
59  clean:
60  	rm -f $(TARGET32).o $(TARGET32) 
\endcode
 *
 * The instructions below describe the processes to make \c Makefile like above.
 *
 *  -# Prepare the environment for making the executable.
 *          
 *      -# Set the directory information. 
 *
 *          In lines 10-14, we set the directories for include files and library file of <b>V</b>.
 *          
 *      -# Set the compiler flags.
 *          
 *          In lines 20-22, we select the compiler, \c g++, set the compiler flags, <tt>-Wall -O2 -m32</tt>,
 *          and set the directory for include files, <tt>-I\$(V_INCDIR)</tt>. 
 *      
 *      -# Set the linker flags.
 *      
 *          In lines 28-31, we set the linker flags to link the target with <b>V</b> library.
 *      
 *      -# Set the source file and the executable file name (lines 38-39).
 *
 *          We set the source file \c myapp.cpp and the executable file name \c myapp.
 *      
 *  -# Create the executable.
 *      
 *      By lines 46-53, \c g++ creates the object file, \c myapp.o, for \c myapp.cpp,
 *      links \c myapp.o with <b>V</b> library, \c libV.a,
 *      and creates the executable file \c myapp.
 *
 *  -# Define the cleaning rule.
 *      
 *      In lines 59-60, we define the rule to remove the object file and the executable. We can remove the object file and the executable by the following:
 *
 *      <tt>\$ make clean</tt>
 * 
 * After making the Makefile, just run \c make and run the executable as the following.
 *      
 * <tt>\$ make</tt>
 *      
 * <tt>\$ ./myapp</tt>
 *
 *
 *
 *
 *
 *
 *
 *
 */







