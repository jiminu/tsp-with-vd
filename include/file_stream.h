#pragma once

#include"city.h"
#include<string>
#include<vector>
#include<map>

using std::string;
using std::vector;
using std::map;
using std::pair;
using distance = float; // == typedef float distance

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
    
    vector<City>& get_cities();
    vector<vector<float>>& get_distnance_matrix();

private :
    // method
    City split_xy(const string& str);
};