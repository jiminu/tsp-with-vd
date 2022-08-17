#include "euler.h"

void Euler::dfs(int u){
    while(1){
        while(!m_pathVector[u].empty() && m_check[m_pathVector[u].top().second]) {
            m_pathVector[u].pop(); // 이미 쓰여진 간선이면 pop
        }
        if(m_pathVector[u].empty()) break;

        auto temp = m_pathVector[u].top(); 
        m_pathVector[u].pop();
        m_check[temp.second] = 1, dfs(temp.first); // y번 간선을 사용했음을 표시하고 알고리즘을 계속한다.
    }

    m_path.push_back(u);
}

Euler::Euler(const vector<pair<int,int>>& edges, const int& n) {
    // vector<vector<int>> adjMatrix;
    // int adjMatrix[100000][100000] = { 0 };
    m_size = n;
    
    // int** adjMatrix = new int*[n];
    // for (int i = 0; i < n; ++i) {
    //     adjMatrix[i] = new int[n];
    // }
    
    m_pathVector.resize(edges.size());
    // adjMatrix.resize(n);
    // for (auto& row : adjMatrix) {
    //     row.resize(n);   
    // }
    int id = 0;
    // for ( int i = 0; i < n; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         adjMatrix[i][j] = 0;
    //     }
    // }
    
    map<pair<int, int>, int> adj;
    for (auto& edge : edges) {
        // if (adjMatrix[edge.first][edge.second] < 0) {
        //     adjMatrix[edge.first][edge.second] = 1;
        //     adjMatrix[edge.second][edge.first] = 1;
            
        // }
        // else {
        //     adjMatrix[edge.first][edge.second]++;
        //     adjMatrix[edge.second][edge.first]++;
        // }
        
            // adjMatrix[edge.first][edge.second]++;
            // adjMatrix[edge.second][edge.first]++;
            adj[{edge.first, edge.second}]++;
    }

    // for(int i = 0; i < n; i++){
    // 	for(int j = i+1; j < n; j++){
    // 		while(adjMatrix[i][j]){
    // 			adjMatrix[i][j]--;
    //             id++;
    // 			m_pathVector[i].push({j, id});
    //             m_pathVector[j].push({i, id}); 
    // 		}
    // 	}
    // }
    for (auto it : adj) {
        id++;
        m_pathVector[it.first.first].push({it.first.second, id});
        m_pathVector[it.first.second].push({it.first.first, id});
    }
    m_check.resize(id + 1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, n);

    dfs(dis(gen));
}

vector<int> Euler::get_path() { return m_path; }