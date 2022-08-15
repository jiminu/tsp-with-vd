#pragma once

#include <vector>

using std::vector;
using std::pair;

class Hamiltonian {
    private:
        vector<int> m_path;
    
        
    public:
        Hamiltonian() {};
        Hamiltonian(const vector<int>& eulerPath);
        vector<int> get_path();
        
    private:
};