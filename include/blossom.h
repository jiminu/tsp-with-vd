#pragma once

#include "Matching.h"
#include <fstream>
#include "Graph.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

class Blossom {
    public:
        pair< Graph, vector<double> > CreateRandomGraph();
        Graph ReadGraph(string filename);
        pair<Graph, vector<double> > ReadWeightedGraph(const int& nodes, vector<pair<pair<int,int>, double>> distanceEdges);
        vector<pair<int, int>> MinimumCostPerfectMatchingExample(const int& nodes, vector<pair<pair<int,int>, double>> distanceEdges);
        void MaximumMatchingExample(string filename);
        // int run(int argc, char* argv[]);
};

