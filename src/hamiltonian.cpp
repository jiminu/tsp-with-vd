#include "hamiltonian.h"
#include "iostream"

#include <map>

using std::map;

Hamiltonian::Hamiltonian(const vector<int>& eulerPath) {
    map<int, int> check;
    
    for (auto node : eulerPath) {
        if (check[node]) {
            continue;
        }
        check[node] = 1;
        m_path.push_back(node);
    }
    
    // m_path.push_back(0);
}

vector<int> Hamiltonian::get_path() { return m_path; }
        