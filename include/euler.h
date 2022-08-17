#pragma once

#include <bits/stdc++.h>
using namespace std;

class Euler {
    private: 
        vector<stack<pair<int, int>>> m_pathVector;
        vector<int> m_check;
        vector<int> m_path;
        int **m_adjMatrix;
        int m_size;

    public:
        Euler() {};
        Euler(const vector<pair<int,int>>& edges, const int& n);
        ~Euler() {
            // for (int i = 0; i < m_size; ++i) {
            //     delete [] m_adjMatrix[i];
            // }
            // delete [] m_adjMatrix;
        }
    
        vector<int> get_path();
        
        
    private: 
        void dfs(int u);
    
};