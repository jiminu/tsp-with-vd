#include "euler.h"

Euler::Euler(int v) {
    vertex = v;
    adj = vector<unordered_map<int, int>>(v + 1);
}

// add edge (u, v) to Euler
void Euler::addEdge(int u, int v) {
    adj[u][v] = 1;
    adj[v][u] = 1;
}

// remove edge (u, v) from the Euler
void Euler::removeEdge(int v, int u) {
    adj[v].erase(u);
    adj[u].erase(v);
}

// function checks if the Euler contains a euler path/circuit or not
// if the Euler is valid, then it calls another function printEuler()
// to print Euler Path or circuit
void Euler::printEulerPathCircuit() {
    int odd = 0;        // number of vertices with odd degree
    int oddVertex = 0;  // it stores vertex with odd degree if it exists

    for (int i = 1; i <= vertex; ++i) {
        if (adj[i].size() % 2 == 1) {
            ++odd;
            oddVertex = i;
        }
    }

    if (odd == 0) {
        // if the number of odd degree vertices is 0
        // the Euler contains a Euler Circuit
        // we can use any vertex as the starting vertex
        cout << "Euler Circuit: ";
        // print euler circuit with '1' as the starting vertex
        printEuler(1);
    } else if (odd == 2) {
        // if the number of odd degree vertices is 0
        // the Euler contains a Euler Path
        // starting vertex should be of odd degree
        cout << "Euler Path: ";
        printEuler(oddVertex);
    } else {
        cout << "Euler Path/Circuit Doesn't Exist" << endl;
    }
}

// the function to print euler path or circuit
void Euler::printEuler(int v) {
    stack<int> cpath;  // current path
    stack<int> epath;  // euler path

    cpath.push(v);  // euler path starts from v

    while (!cpath.empty()) {
        int u = cpath.top();

        if (adj[u].size() == 0) {
            // if all edges from u are visited
            // pop u and push it to euler path
            epath.push(u);
            cpath.pop();
        } else {
            // if all edges from u are not visited
            // select any edge (u, v)
            // push v to current path and remove edge (u, v) from the Euler
            cpath.push(adj[u].begin()->first);
            removeEdge(u, adj[u].begin()->first);
        }
    }

    while (!epath.empty()) {
        cout << " " << epath.top() << " ";
        epath.pop();
    }
}