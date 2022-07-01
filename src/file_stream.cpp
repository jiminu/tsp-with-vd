#include"file_stream.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<math.h>

#include <random>
#include<algorithm>
#include <VoronoiDiagram2DC.h>

FileStream::FileStream() {
}

void FileStream::read(const string& fileName) {
    std::ifstream file(fileName);
    int run = 0;
    string line;
    
    if(file) {
        while(std::getline(file, line)) {
            if (file.eof()) break;
            if (run < 6) {
                ++run;
                continue;
            }
            m_citiesVector.push_back(split_xy(line));
        }
        file.close();
    }
    else {
        std::cout << "file open fail" << std::endl;
        exit(0);
    }
}

// void FileStream::read(const string& fileName) {
//     std:: ifstream file(fileName);
//     bool firstRun = true;
//     string line;
    
//     if(file) {
//         while(std::getline(file, line)) {
//             if (file.eof()) break;
//             if (firstRun) {
//                 firstRun = false;
//                 continue;
//             }
//             m_citiesVector.push_back(split_xy(line));
//         }
//         file.close();
//     }
//     else {
//         std::cout << "file open fail" << std::endl;
//         exit(0);
//     }
// }

void FileStream::read_distance_matrix(const string& fileName) {
    std::ifstream file(fileName);
    int run = 0;
    string line;
    
    if(file) {
        while(std::getline(file, line)) {
            // if (file.eof()) break;
            if (run < 1) {
                ++run;
                continue;
            }
            std::istringstream ss(line);
            string stringBuffer;
            vector<float> tempVector;
            while (getline(ss, stringBuffer, ',')) {
                tempVector.push_back(std::stof(stringBuffer));
            }
            m_distanceMatrix.push_back(tempVector);
        }
        file.close();
    }
    else {
        std::cout << "file open fail" << std::endl;
        exit(0);
    }
}

void FileStream::write(const string& fileName, const pair<float, vector<int>>& bestSolution, const vector<float>& info) {
    std::ofstream fout(fileName);

    fout << "fitness : " << bestSolution.first << "\n";
    fout << "selection pressure : " << info[0] << "\n";
    fout << "crossover parameter : " << info[1] << "\n";
    fout << "mutation parameter : " << info[2] << "\n";
    fout << "population : " << info[3] << "\n";
    fout << "elite proportion : " << info[4] << "\n";
    fout << "generation : " << info[5] << "\n";
    fout << "time : " << info[6] << "\n";
    for (auto it = bestSolution.second.begin(); it != bestSolution.second.end(); it++) {
        fout << *it << "\n";
    }
    
    fout.close();
}

void FileStream::write_to_edge(const string& fileName, const list<VEdge2D*>& edge) {
    std::ofstream fout("./../data/answer_voronoi_edge.txt");
    for (const auto& it : edge) {
        fout << it->getStartVertex()->getLocation().getX() << "," << it->getStartVertex()->getLocation().getY() << "," <<
        it->getEndVertex()->getLocation().getX() << "," << it->getEndVertex()->getLocation().getY() << "\n";
    }
    fout.close();
}

void FileStream::write_to_face(const string& fileName, const list<VFace2D*>& face) {
    std::ofstream fout("./../data/answer_voronoi_face.txt");
    for (const auto& it : face) {
        fout << it->getGenerator()->getDisk().getX() << "," << it->getGenerator()->getDisk().getY() << "," <<
        it->getGenerator()->getDisk().getX() << "," << it->getGenerator()->getDisk().getY() << "\n";
    }
    fout.close();
}

void FileStream::write_to_vertices(const string& fileName, const list<VVertex2D*>& vertices) {
    std::ofstream fout("./../data/answer_voronoi_vertices.txt");
    for (const auto& it : vertices) {
        fout << it->getLocation().getX() << "," << it->getLocation().getY() << "," <<
        it->getLocation().getX() << "," << it->getLocation().getY() << "\n";
    }
    fout.close();
}

// City FileStream::split_xy(const string& str) {
//     City city;
    
//     city.x = std::stof(str.substr(4, 3));
//     city.y = std::stof(str.substr(8, 3));
    
//     return city;
// }

City FileStream::split_xy(const string& str)  {
    std::istringstream split_str(str); 
    string buffer;
    City city;
    int i = 0;
    while(std::getline(split_str, buffer, ' ')) {
        if(i == 1) {
            city.x = stof(buffer);
        }
        else if(i == 2) {
            city.y = stof(buffer);
        }
        ++i;
    }
    return city;
}

vector<City>& FileStream::get_cities() {
    return m_citiesVector;
}

vector<vector<float>>& FileStream::get_distnance_matrix() {
    return m_distanceMatrix;
}