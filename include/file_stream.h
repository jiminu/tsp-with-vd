#pragma once

#include"city.h"
#include<string>
#include<vector>
#include<map>
#include<list>

#include<VEdge2D.h>
#include<VFace2D.h>

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::list;
using distance = float; // == typedef float distance

using namespace V::GeometryTier;

class FileStream {
private :
    // member variable
    vector<City> m_citiesVector;
    vector<vector<float>> m_distanceMatrix;

public :
    // method
    FileStream();
    ~FileStream() {};
    
    void read(const string& fileName);
    void read_distance_matrix(const string& fileName);
    void write(const string& fileName, const pair<float, vector<int>>& bestSolution, const vector<float>& info);
    
    // void write_to_vertex(const string& fileName, const vector<Vertex_s*>& vertices);
    void write_to_edge(const string& fileName, const list<VEdge2D*>& edge);
    void write_to_face(const string& fileName, const list<VFace2D*>& face);
    
    vector<City>& get_cities();
    vector<vector<float>>& get_distnance_matrix();

private :
    // method
    City split_xy(const string& str);
    bool check_same_chain(VFace2D* currFace, VFace2D* targetFace, vector<list<VFace2D*>>& connectedChain);
    bool check_target_face_state(VFace2D* targetFace, map<VFace2D*, int>& connectedFaces);
    void connect_chain(VFace2D* currFace,
                       VFace2D* targetFace,
                       multimap<int, VFace2D*>& chainCountEdges,
                       vector<list<VFace2D*>>& connectedChain,
                       map<VFace2D*, int>& connectedFaces
    );
};