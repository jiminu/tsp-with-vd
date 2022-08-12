#include <bits/stdc++.h>
using namespace std;

class Eul {
    private: 
        vector<stack<pair<int, int>>> pathVector;
        vector<int> check;
        vector<int> result;

    public:
        void dfs(int u);

        vector<int> run(const vector<pair<int,int>>& edges, const int& n);
};