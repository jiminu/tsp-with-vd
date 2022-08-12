#include "hamiltonian.h"
#include "iostream"

#include <map>

using std::map;

vector<int> Hamiltonian::run(const vector<int>& eulerPath) {
    vector<int> result;
    map<int, int> check;
    
    for (auto node : eulerPath) {
        if (check[node]) {
            continue;
        }
        check[node] = 1;
        result.push_back(node);
    }
    
    for (auto a : result) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
    std::cout << result.size() << std::endl;
    
    
    
    return result;
}
