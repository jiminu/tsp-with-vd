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

// void FileStream::read(const string& fileName) {
//     std::ifstream file(fileName);
//     int run = 0;
//     string line;
    
//     if(file) {
//         while(std::getline(file, line)) {
//             if (file.eof()) break;
//             if (run < 6) {
//                 ++run;
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

void FileStream::read(const string& fileName) {
    std:: ifstream file(fileName);
    bool firstRun = true;
    string line;
    
    if(file) {
        while(std::getline(file, line)) {
            if (file.eof()) break;
            if (firstRun) {
                firstRun = false;
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
    int gg = 0;
    std::ofstream fout("./../data/answer_voronoi_face.txt");
    int count = 0;
    bool status = false;
    multimap<int, VFace2D*> chainCountEdges;
    vector<list<VFace2D*>> connectedChain;
    map<VFace2D*, int> connectedFaces;
    
    vector<VFace2D*> faceVector(face.begin(), face.end());
    // std::random_device rd;
    // std::mt19937 g(rd());
    // std::shuffle(faceVector.begin(), faceVector.end(), g);
    
    
    for (const auto& face : faceVector) {
        chainCountEdges.insert({0, face});
    }
    
    std::cout << "0 : " << chainCountEdges.count(0) << std::endl;
    std::cout << "1 : " << chainCountEdges.count(1) << std::endl;
    std::cout << "2 : " << chainCountEdges.count(2) << std::endl;
    std::cout << "3 : " << chainCountEdges.count(3) << "\n" << std::endl;

    auto currFaceIter = chainCountEdges.begin();
    auto currFace = chainCountEdges.begin()->second;
    while (chainCountEdges.count(0) > 0 || chainCountEdges.count(1) > 0) {        
        gg++;
        list<VEdge2D*> boundaryEdges;
        currFace->getBoundaryVEdges(boundaryEdges);
        map<float, VFace2D*> boundaryFaces;
        for (const auto& boundaryEdge : boundaryEdges) {
            if (boundaryEdge->getRightFace()->getGenerator() == 0 ||
                boundaryEdge->getLeftFace()->getGenerator() == 0 ||
                boundaryEdge->getRightFace()->getGenerator()->getID() == -1 ||
                boundaryEdge->getLeftFace()->getGenerator()->getID() == -1) {
                continue;
            }

            float dist = boundaryEdge->getLeftFace()->getGenerator()->getDisk().getCenterPt()
                         .distance(boundaryEdge->getRightFace()->getGenerator()->getDisk().getCenterPt());

            if (boundaryEdge->getLeftFace() == currFace) {
                boundaryFaces.insert({dist, boundaryEdge->getRightFace()});
            }

            else {
                boundaryFaces.insert({dist, boundaryEdge->getLeftFace()});
            }
        }
        
        for (auto target = boundaryFaces.begin(); target != boundaryFaces.end(); ++target) {
            VFace2D* targetFace = target->second;

            if (check_same_chain(currFace, targetFace, connectedChain) || check_target_face_state(targetFace, connectedFaces)) continue;
            
            connect_chain(currFace, targetFace, chainCountEdges, connectedChain, connectedFaces);
            fout << currFace->getGenerator()->getDisk().getX() << "," << currFace->getGenerator()->getDisk().getY() << ","
                 << targetFace->getGenerator()->getDisk().getX() << "," << targetFace->getGenerator()->getDisk().getY() << "\n";
            // std::cout << currFace->getGenerator()->getDisk().getX() << "," << currFace->getGenerator()->getDisk().getY() << ","
            //      << targetFace->getGenerator()->getDisk().getX() << "," << targetFace->getGenerator()->getDisk().getY() << "\n";
            ++count;
            status = true;
            currFaceIter = chainCountEdges.begin();
            currFace = currFaceIter->second;
            break;
        }
        
        if (!status) {
            map<float, VFace2D*> tempMap;
            auto candidate = chainCountEdges.lower_bound(1);
            candidate++;
            for (; candidate != chainCountEdges.upper_bound(1); ++candidate) {
                float dist = currFace->getGenerator()->getDisk().getCenterPt().distance(candidate->second->getGenerator()->getDisk().getCenterPt());
                tempMap.insert({dist, candidate->second});
            }
            auto tempFace = tempMap.begin();
            while (check_same_chain(currFace, tempFace->second, connectedChain) || check_target_face_state(tempFace->second, connectedFaces)) {
                tempFace++;
            }
            connect_chain(currFace, tempFace->second, chainCountEdges, connectedChain, connectedFaces);
            fout << currFace->getGenerator()->getDisk().getX() << "," << currFace->getGenerator()->getDisk().getY() << ","
                 << tempFace->second->getGenerator()->getDisk().getX() << "," << tempFace->second->getGenerator()->getDisk().getY() << "\n";
            currFaceIter = chainCountEdges.begin();
            currFace = currFaceIter->second;
        }
        status = false;
        if (chainCountEdges.count(2) > 1 && connectedChain.size() == 1) break;
        
        
        // if (currFaceIter == chainCountEdges.end()) {
        //     currFaceIter = chainCountEdges.begin();
        //     map<float, VFace2D*> tempMap;
            
        //     for (auto candidate = chainCountEdges.lower_bound(1); candidate != chainCountEdges.upper_bound(1); ++candidate) {
        //         float dist = currFace->getGenerator()->getDisk().getCenterPt()
        //                      .distance(candidate->second->getGenerator()->getDisk().getCenterPt());

        //         tempMap.insert({dist, candidate->second});
        //     }
            
        //     auto tempFace = tempMap.begin();
        //     while (check_same_chain(currFace, tempFace->second, connectedChain) || check_target_face_state(tempFace->second, connectedFaces)) {
        //         tempFace++;
        //     }
        //     connect_chain(currFace, tempFace->second, chainCountEdges, connectedChain, connectedFaces);
        //     fout << currFace->getGenerator()->getDisk().getX() << "," << currFace->getGenerator()->getDisk().getY() << ","
        //          << tempFace->second->getGenerator()->getDisk().getX() << "," << tempFace->second->getGenerator()->getDisk().getY() << "\n";
            
        // }
        // else if (!status) {
        //     currFaceIter++;
        //     currFace = currFaceIter->second;
        // }
        // status = false;
    }
    
    // ===============================================================================================================
    // while (chainCountEdges.count(0) > 0) {
    //     VFace2D* currFace = chainCountEdges.lower_bound(0)->second;

    //     list<VEdge2D*> boundaryEdges;
    //     currFace->getBoundaryVEdges(boundaryEdges);
    //     map<float, VFace2D*> boundaryFaces;

    //     for (const auto& boundaryEdge : boundaryEdges) {
    //         if (boundaryEdge->getRightFace()->getGenerator() == 0 ||
    //             boundaryEdge->getLeftFace()->getGenerator() == 0) {
    //             continue;
    //         }

    //         float dist = boundaryEdge->getLeftFace()->getGenerator()->getDisk().getCenterPt()
    //                       .distance(boundaryEdge->getRightFace()->getGenerator()->getDisk().getCenterPt());

    //         if (boundaryEdge->getLeftFace() == currFace) {
    //             boundaryFaces.insert({dist, boundaryEdge->getRightFace()});
    //         } 
            
    //         else {
    //             boundaryFaces.insert({dist, boundaryEdge->getLeftFace()});
    //         }
    //     }
        
    //     for (auto target = boundaryFaces.begin(); target != boundaryFaces.end(); ++target) {
    //         VFace2D* targetFace = target->second;

    //         if (check_same_chain(currFace, targetFace, connectedChain)) continue;
    //         if (check_target_face_state(targetFace, connectedFaces)) continue;
            
    //         connect_chain(currFace, targetFace, chainCountEdges, connectedChain, connectedFaces);
    //         fout << currFace->getGenerator()->getDisk().getX() << "," << currFace->getGenerator()->getDisk().getY() << ","
    //              << targetFace->getGenerator()->getDisk().getX() << "," << targetFace->getGenerator()->getDisk().getY() << "\n";
    //         ++count;
    //         break;
    //     }
    // }
    
    // std::cout << "0 connected done : " << count << std::endl;
    // std::cout << "0 : " << chainCountEdges.count(0) << std::endl;
    // std::cout << "1 : " << chainCountEdges.count(1) << std::endl;
    // std::cout << "2 : " << chainCountEdges.count(2) << std::endl;
    // std::cout << "3 : " << chainCountEdges.count(3) << "\n" << std::endl;
        
    // while (chainCountEdges.count(1) > 0) {
    //     VFace2D* currFace = chainCountEdges.lower_bound(1)->second;

    //     list<VEdge2D*> boundaryEdges;
    //     currFace->getBoundaryVEdges(boundaryEdges);
    //     map<float, VFace2D*> boundaryFaces;

    //     for (const auto& boundaryEdge : boundaryEdges) {
    //         if (boundaryEdge->getRightFace()->getGenerator() == 0 ||
    //             boundaryEdge->getLeftFace()->getGenerator() == 0) {
    //             continue;
    //         }

    //         float dist = boundaryEdge->getLeftFace()->getGenerator()->getDisk().getCenterPt()
    //                       .distance(boundaryEdge->getRightFace()->getGenerator()->getDisk().getCenterPt());

    //         if (boundaryEdge->getLeftFace() == currFace) {
    //             boundaryFaces.insert({dist, boundaryEdge->getRightFace()});
    //         }

    //         else {
    //             boundaryFaces.insert({dist, boundaryEdge->getLeftFace()});
    //         }
    //     }

    //     for (auto target = boundaryFaces.begin(); target != boundaryFaces.end(); ++target) {
    //         VFace2D* targetFace = target->second;
            
    //         if (targetFace == boundaryFaces.cbegin()->second && (check_same_chain(currFace, targetFace, connectedChain) || check_target_face_state(targetFace, connectedFaces))) {
    //             map<float, VFace2D*> temp;
    //             for (auto candidate = ++chainCountEdges.lower_bound(1); candidate != chainCountEdges.upper_bound(1); ++candidate) {
    //                 float dist = currFace->getGenerator()->getDisk().getCenterPt()
    //                              .distance(candidate->second->getGenerator()->getDisk().getCenterPt());
    //                 temp.insert({dist, candidate->second});
    //             }
                
    //             auto tempFace = temp.begin();
    //             while (check_same_chain(currFace, tempFace->second, connectedChain) || check_target_face_state(tempFace->second, connectedFaces)) {
    //                 tempFace++;
    //             }
                
                    
    //             connect_chain(currFace, tempFace->second, chainCountEdges, connectedChain, connectedFaces);
    //             fout << currFace->getGenerator()->getDisk().getX() << "," << currFace->getGenerator()->getDisk().getY() << ","
    //                  << tempFace->second->getGenerator()->getDisk().getX() << "," << tempFace->second->getGenerator()->getDisk().getY() << "\n";
    //             ++count;
    //             break;
    //         }
    //         else if (check_same_chain(currFace, targetFace, connectedChain) || check_target_face_state(targetFace, connectedFaces)) {
    //             continue;
    //         }
    //         connect_chain(currFace, targetFace, chainCountEdges, connectedChain, connectedFaces);
    //         fout << currFace->getGenerator()->getDisk().getX() << "," << currFace->getGenerator()->getDisk().getY() << ","
    //              << targetFace->getGenerator()->getDisk().getX() << "," << targetFace->getGenerator()->getDisk().getY() << "\n";
    //         ++count;
    //         break;
    //     }
    // }
    // ===============================
    fout.close();
    
    // std::ofstream fout2("./../data/checkPoint.txt");
    // for (auto it = chainCountEdges.lower_bound(1); it != chainCountEdges.upper_bound(1); ++it) {
    //     fout2 << it->second->getGenerator()->getDisk().getX() << "," << it->second->getGenerator()->getDisk().getY() << "\n";
    // }
    // fout2.close();
    
    std::cout << "1 connected done : " << count << std::endl;
    std::cout << "0 : " << chainCountEdges.count(0) << std::endl;
    std::cout << "1 : " << chainCountEdges.count(1) << std::endl;
    std::cout << "2 : " << chainCountEdges.count(2) << std::endl;
    std::cout << "3 : " << chainCountEdges.count(3) << "\n" << std::endl;
}

bool FileStream::check_same_chain(VFace2D* currFace, VFace2D* targetFace, vector<list<VFace2D*>>& connectedChain) {
    for (const auto& currentChain : connectedChain) {
        if ((currentChain.front() == currFace && currentChain.back() == targetFace) ||
            (currentChain.front() == targetFace && currentChain.back() == currFace)) {
                return true;
        }
    }
    return false;
}

bool FileStream::check_target_face_state(VFace2D* targetFace, map<VFace2D*, int>& connectedFaces) {
    if (connectedFaces[targetFace] == 2) {
        return true;
    }
    return false;
}

void FileStream::connect_chain(VFace2D* currFace,
                               VFace2D* targetFace,
                               multimap<int, VFace2D*>& chainCountEdges,
                               vector<list<VFace2D*>>& connectedChain,
                               map<VFace2D*, int>& connectedFaces) {
    int currentConnectedNode = connectedFaces[currFace];
    int targetConnectedNode = connectedFaces[targetFace];
    connectedFaces[currFace] = currentConnectedNode + 1;
    connectedFaces[targetFace] = targetConnectedNode + 1;
    
    for (auto it = chainCountEdges.lower_bound(currentConnectedNode); it != chainCountEdges.upper_bound(currentConnectedNode); ++it) {
        if (currFace == it->second) {
            chainCountEdges.erase(it);
            break;
        }
    }
    for (auto it = chainCountEdges.lower_bound(targetConnectedNode); it != chainCountEdges.upper_bound(targetConnectedNode); ++it) {
        if (targetFace == it->second) {
            chainCountEdges.erase(it);
            break;
        }
    }
    
    chainCountEdges.insert({currentConnectedNode+1, currFace});
    chainCountEdges.insert({targetConnectedNode+1, targetFace});
    
    string currPos = "none";
    string targetPos = "none";
    
    int startPos;
    int endPos;
    
    for (int i = 0; i < connectedChain.size(); ++i) {
        if (currFace == connectedChain[i].front()) {
            currPos = "front";
            startPos = i;
        }
        else if(currFace == connectedChain[i].back()) {
            currPos = "back";
            startPos = i;
        }
        if (targetFace == connectedChain[i].front()) {
            targetPos = "front";
            endPos = i;
        } 
        else if(targetFace == connectedChain[i].back()) {
            targetPos = "back";
            endPos = i;
        }
    }
    
    if (currPos == "none" && targetPos == "none") {
        connectedChain.push_back({currFace, targetFace});
    }
    else if (currPos == "none") {
        if (targetPos == "front") {
            connectedChain[endPos].push_front(currFace);
        }
        else {
            connectedChain[endPos].push_back(currFace);
        }
    }
    else if (targetPos == "none") {
        if (currPos == "front") {
            connectedChain[startPos].push_front(targetFace);            
        }
        else {
            connectedChain[startPos].push_back(targetFace);            
        }
    }
    else {
        if (currPos == "front" && targetPos == "back") {
            connectedChain[endPos].insert(connectedChain[endPos].end(), connectedChain[startPos].begin(),connectedChain[startPos].end());
            connectedChain.erase(connectedChain.begin() + startPos);
        }
        else if (currPos == "front" && targetPos == "front") {
            connectedChain[startPos].reverse();
            connectedChain[startPos].insert(connectedChain[startPos].end(), connectedChain[endPos].begin(),connectedChain[endPos].end());
            connectedChain.erase(connectedChain.begin() + endPos);
        }
        else if (currPos == "back" && targetPos == "front") {
            connectedChain[startPos].insert(connectedChain[startPos].end(), connectedChain[endPos].begin(),connectedChain[endPos].end());
            connectedChain.erase(connectedChain.begin() + endPos);
        }
        else if (currPos == "back" && targetPos == "back") {
            connectedChain[startPos].reverse();
            connectedChain[endPos].insert(connectedChain[endPos].end(), connectedChain[startPos].begin(),connectedChain[startPos].end());
            connectedChain.erase(connectedChain.begin() + startPos);
        }
    }
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
    while(std::getline(split_str, buffer, '\t')) {
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