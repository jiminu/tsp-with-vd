#pragma once

#include<string>
#include<vector>
#include<map>
#include<ctime>
#include "city.h"
#include <VoronoiDiagram2DC.h>
#include <VoronoiDiagramCIC.h>

using std::vector;
using std::map;
using std::multimap;
using std::pair;
using std::string;

using namespace V::GeometryTier;

class HeuristicAlgorithm {
    private:
        
        const float m_selectionPressure  = 3;
        const float m_eliteProportion    = 0.2;
        const float m_crossoverParameter = 0.7;
        const float m_mutationParameter  = 0.2;
        const int m_population           = 100;
        const int m_generation           = 10000000;
        const string m_mutation          = "inversion";
        
        vector<vector<float>> m_distanceMatrix;

        string m_tspFile  = "./../data/tsp_data.txt";
        string m_distanceMatrixFile  = "./../data/dist.txt";
        
        string m_savePath = "./../data/";
        string m_saveFile = m_mutation + ".txt";
    
        // VoronoiDiagram2DC m_VD;
        VoronoiDiagramCIC m_VD;
        vector<City> m_cities;
        pair<float, vector<int>> m_bestSolution = {0, {}};
        int m_currGeneration = 0;

        clock_t start, end;
        float result;
    public:
    
        HeuristicAlgorithm();  
        ~HeuristicAlgorithm();  
        
        
    private:
        vector<pair<float, vector<int>>> initialize_chromosome(const int& population);
        vector<pair<float, vector<int>>> selection(const vector<pair<float, vector<int>>>& chromosomes);
        vector<pair<float, vector<int>>> evaluation(const vector<vector<int>>& populations);
        vector<pair<float, vector<int>>> crossover(vector<pair<float, vector<int>>>& selectionPopulations);
        void mutation(pair<float, vector<int>>& offspring);
        
        void generate_vd();
        vector<pair<float, vector<int>>> initialize_chromosome_with_VD(const int& population);
        
        
        void order_crossover(vector<pair<float, vector<int>>>& selectionPopulations);
        
        void insertion_mutation(pair<float, vector<int>>& crossoverPopulations);
        void inversion_mutation(pair<float, vector<int>>& crossoverPopulations);
        void displacement_mutation(pair<float, vector<int>>& crossoverPopulations);

        void check_same_value(vector<int>& edge, const int& value);
        void erase_value_from_edge(map<int, vector<int>>& edge, const int& value);
        vector<pair<float, vector<int>>*> select_parents(vector<pair<float, vector<int>>>& selectionPopulations);
        float evaluate_function(const vector<int>& population);
        void save_best_solution(const vector<float>& info);

        pair<float, vector<int>> find_best_fitness(const vector<pair<float, vector<int>>>& populations);

        void generate_cities();
        void generate_distance_matrix();
        int generate_random_int(const int& min, const int& max);
        float generate_random_float(const float& min, const float& max);
};