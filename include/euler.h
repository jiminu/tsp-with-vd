#pragma once

#include <bits/stdc++.h>
using namespace std;

class Euler {
    private: 
        vector<stack<pair<int, int>>> m_pathVector;
        vector<int> m_check;
        vector<int> m_path;

    public:
        Euler() {};
        Euler(const vector<pair<int,int>>& edges, const int& n);
    
        vector<int> get_path();
        
        
    private: 
        void dfs(int u);
    
};