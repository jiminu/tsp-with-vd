#pragma once

#include<iostream>
#include<vector>
#include<stack>
#include<unordered_map>
using namespace std;

class Euler{
    
    int vertex; // number of vertices
    // adjacency list
    // we are using unordered_map to speed up the deletion of edges
    vector<unordered_map<int,int>> adj;
    
    public:
        
        // constructor to initialize Euler
        // set number of vertices to v
        // initialize adjacency map for v nodes
        Euler(int v);
        
        // add edge (u, v) to Euler
        void addEdge(int u, int v);
        
        // remove edge (u, v) from the Euler
        void removeEdge(int v,int u);
        
        // function checks if the Euler contains a euler path/circuit or not
        // if the Euler is valid, then it calls another function printEuler()
        // to print Euler Path or circuit
        void printEulerPathCircuit();
        
        // the function to print euler path or circuit
        void printEuler(int v);
        
};