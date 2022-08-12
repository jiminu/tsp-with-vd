#include "eul.h"

void Eul::dfs(int u){
    while(1){
        while(!pathVector[u].empty() && check[pathVector[u].top().second]) {
            pathVector[u].pop(); // 이미 쓰여진 간선이면 pop
        }
        if(pathVector[u].empty()) break;

        auto temp = pathVector[u].top(); 
        pathVector[u].pop();
        check[temp.second] = 1, dfs(temp.first); // y번 간선을 사용했음을 표시하고 알고리즘을 계속한다.
    }

    result.push_back(u);
    // std::cout << u << " ";
}

vector<int> Eul::run(const vector<pair<int,int>>& edges, const int& n) {
    vector<vector<int>> adjMatrix;
    
    pathVector.resize(edges.size());
    adjMatrix.resize(n);
    for (auto& row : adjMatrix) {
        row.resize(n);   
    }
    int id;
    
    for (auto& edge : edges) {
        adjMatrix[edge.first][edge.second]++;
        adjMatrix[edge.second][edge.first]++;
    }

    for(int i = 0; i < n; i++){
    	for(int j = i+1; j < n; j++){
    		while(adjMatrix[i][j]){
    			adjMatrix[i][j]--;
                id++;
    			pathVector[i].push({j, id}), pathVector[j].push({i, id}); // 간선에 번호를 부여한다.
    		}
    	}
    }
    check.resize(id + 1);

    // for(int i = 0; i < n; i++) {
    //     if(pathVector[i].size() % 2) {
    //         return !printf("-1");
    //     }
    // }

    dfs(0);
    
    return result;
}