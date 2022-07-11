#include "heuristic_algorithm.h"
#include <time.h>

#include <algorithm>
#include <random>
#include <iostream>

#include "file_stream.h"
#include "QuasiTriangulation2D.h"

HeuristicAlgorithm::HeuristicAlgorithm() {
    start = clock();
    generate_cities();
    // generate_distance_matrix();
    generate_vd();

    float start = clock();
    
    // vector<pair<float, vector<int>>> populations = initialize_chromosome_with_VD(m_population);
    // // vector<pair<float, vector<int>>> populations = initialize_chromosome(m_population);
    // float end = clock();
    // vector<float> info = {m_selectionPressure, 
    //                       m_crossoverParameter, 
    //                       m_mutationParameter, 
    //                       float(m_population), 
    //                       float(m_generation),
    //                       m_eliteProportion, 
    //                       result / CLOCKS_PER_SEC};
    // save_best_solution(info);
    // std::cout << "initialize chromosome time : " << (end - start) / CLOCKS_PER_SEC << "s" << std::endl;

    // for (auto it : populations) {
    //     std::cout << it.first << std::endl;
    // }

    // for (int i = 0; i < m_generation; ++i) {
    //     // populations = selection(populations);
    //     // populations = crossover(populations);
        
    //     m_currGeneration = i;
    //     // if (i % 100 == 0) {
    //         std::cout << "Generation " << i << " : " << find_best_fitness(populations).first << std::endl;
    //     // }
    // }
    // end = clock();
    // result = end - start;
    // info = {m_selectionPressure, 
    //                       m_crossoverParameter, 
    //                       m_mutationParameter, 
    //                       float(m_population), 
    //                       float(m_generation),
    //                       m_eliteProportion, 
    //                       result / CLOCKS_PER_SEC};
    
    // std::cout << "fitness : " << m_bestSolution.first << std::endl;
    // for (const auto& it : m_bestSolution.second) {
    //     std::cout << it << " -> ";
    // }
    // std::cout << m_bestSolution.second[0] << std::endl;
    // std::cout << "runtime : " << result / CLOCKS_PER_SEC << "s" << std::endl;
    // std::cout << "best solution : " << m_bestSolution.first << std::endl;
    // save_best_solution(info);
}

HeuristicAlgorithm::~HeuristicAlgorithm() {
    
}

void HeuristicAlgorithm::generate_vd() {
    list<rg_Circle2D> circles;
    
    for (const auto& city : m_cities) {
        rg_Circle2D circle(city.x, city.y, 0);
        
        circles.push_back(circle);
    }
    std::ofstream fout("./../data/delaunay.txt");
    float start = clock();
    m_VD.constructVoronoiDiagram(circles);
    // m_VD.constructVoronoiDiagramCIC_noContainerInInput(circles);
    QuasiTriangulation2D QT;
    QT.construct(m_VD);
    m_BU.construct(QT);
    
    multimap<double, EdgeBU2D> tempMap;
    
    rg_dList<EdgeBU2D> tempList;
    tempList = m_BU.getEdges();
    
    tempList.reset4Loop();
    while ( tempList.setNext4Loop() ) {
        EdgeBU2D currEdge = tempList.getEntity();
        
        if( currEdge.isVirtual() ) {
            continue;
        }
        fout << currEdge.getStartVertex()->getCircle().getX() << "," << currEdge.getStartVertex()->getCircle().getY() << "," 
            << currEdge.getEndVertex()->getCircle().getX() << "," << currEdge.getEndVertex()->getCircle().getY() << "\n";
            
        tempMap.insert({currEdge.getStartVertex()->getCoord().distance(currEdge.getEndVertex()->getCoord()), currEdge});
    }
    fout.close();
    
    generate_mst(tempMap);
    float end = clock();
    
    std::cout << "generate voronoi diagram time : " << (end - start) / CLOCKS_PER_SEC << "s" << std::endl;
    
    string face_path = "./../data/answer_voronoi_faces.txt";
    string edge_path = "./../data/answer_voronoi_edges.txt";
    string vertex_path = "./../data/answer_voronoi_vertices.txt";
    
    list<VFace2D*> faces; 
    list<VEdge2D*> edges; 
    list<VVertex2D*> vertices; 

    m_VD.getVoronoiFaces(faces);
    m_VD.getVoronoiEdges(edges);
    m_VD.getVoronoiVertices(vertices);
    faces.pop_front();
    
    save_face(face_path, faces);
    save_edge(edge_path, edges);
    save_vertex(vertex_path, vertices);
}

int HeuristicAlgorithm::find_parents(vector<int>& set, const int id) {
    if (set[id] == id) return id;
    
    return set[id] = find_parents(set, set[id]);
}
void HeuristicAlgorithm::union_parents(vector<int>& set, int a, int b) {
    a = find_parents(set, a);
    b = find_parents(set, b);
    
    if (a < b) {
        set[b] = a;
    }
    else {
        set[a] = b;
    }
}

void HeuristicAlgorithm::generate_mst(const multimap<double, EdgeBU2D>& distanceMap) {
    std::ofstream fout("./../data/result1.txt");
    std::ofstream fout2("./../data/odd_face.txt");
    std::ofstream fout3("./../data/odd_edge.txt");
    std::ofstream fout4("./../data/odd_mst.txt");
    
    int count = 0;
    
    map<VertexBU2D*, int> connectedFaces;
    vector<int> set;
    list<EdgeBU2D> resultEdges;
    
    for (int i = 0; i < m_cities.size(); ++i) {
        set.push_back(i);
    }
    
    for (auto& edge : distanceMap) {
        int startVertexID = edge.second.getStartVertex()->getID();
        int endVertexID = edge.second.getEndVertex()->getID();
        
        if ( find_parents(set, startVertexID) != find_parents(set, endVertexID)) {
            union_parents(set, startVertexID, endVertexID);
            
            connectedFaces.insert({edge.second.getStartVertex(), connectedFaces[edge.second.getStartVertex()]++});
            connectedFaces.insert({edge.second.getEndVertex(), connectedFaces[edge.second.getEndVertex()]++});
            
            resultEdges.push_back(edge.second);
            fout << edge.second.getStartVertex()->getCircle().getX() << "," << edge.second.getStartVertex()->getCircle().getY() << "," 
                 << edge.second.getEndVertex()->getCircle().getX()  << "," << edge.second.getEndVertex()->getCircle().getY()  << "\n";
            if (resultEdges.size() == m_cities.size() - 1) break;
            continue;
        }
    }
    fout.close();
    
    list<rg_Circle2D> oddFaces;
    
    for (auto it : connectedFaces) {
        if (it.second % 2 != 0) {
             fout2 << it.first->getCircle().getX() << "," << it.first->getCircle().getY() << "\n";
             oddFaces.push_back(it.first->getCircle());
        }
    }
    fout2.close();
    
    // odd degree vertices calculate..
    
    map<VertexBU2D*, int> connectedFaces2;
    vector<int> set2;
    list<EdgeBU2D> resultEdges2;
    
    for (int i = 0; i < m_cities.size(); ++i) {
        set2.push_back(i);
    }
    
    for (auto& edge : distanceMap) {
        int startVertexID = edge.second.getStartVertex()->getID();
        int endVertexID = edge.second.getEndVertex()->getID();
        
        if ( find_parents(set2, startVertexID) != find_parents(set2, endVertexID)) {
            union_parents(set2, startVertexID, endVertexID);
            
            connectedFaces2.insert({edge.second.getStartVertex(), connectedFaces2[edge.second.getStartVertex()]++});
            connectedFaces2.insert({edge.second.getEndVertex(), connectedFaces2[edge.second.getEndVertex()]++});
            
            resultEdges2.push_back(edge.second);
            fout4 << edge.second.getStartVertex()->getCircle().getX() << "," << edge.second.getStartVertex()->getCircle().getY() << "," 
                 << edge.second.getEndVertex()->getCircle().getX()  << "," << edge.second.getEndVertex()->getCircle().getY()  << "\n";
            if (resultEdges2.size() == oddFaces.size() - 1) break;
            continue;
        }
    }
    fout4.close();
    
    // list<VEdge2D*> oddEdges;
    // VoronoiDiagram2DC oddVD;
    // oddVD.constructVoronoiDiagram(oddFaces);
    // list<VFace2D*> faces; 
    // oddVD.getVoronoiFaces(faces);
    
    // map<VFace2D*, int> degFaces;
    
    // for (auto& face : faces) {
    //     degFaces[face] = 0;
    // }
    
    // for (auto& face : degFaces) {
    //     if (face.second != 0) continue;
    //     list<VEdge2D*> boundaryEdges;
    //     face.first->getBoundaryVEdges(boundaryEdges);
    //     map<float, VFace2D*> boundaryFaces;
    //     for (const auto& boundaryEdge : boundaryEdges) {
    //         if (boundaryEdge->getRightFace()->getGenerator() == 0 ||
    //             boundaryEdge->getLeftFace()->getGenerator() == 0 ||
    //             boundaryEdge->getRightFace()->getGenerator()->getID() == -1 ||
    //             boundaryEdge->getLeftFace()->getGenerator()->getID() == -1) {
    //             continue;
    //         }

    //         float dist = boundaryEdge->getLeftFace()->getGenerator()->getDisk().getCenterPt()
    //                      .distance(boundaryEdge->getRightFace()->getGenerator()->getDisk().getCenterPt());

    //         if (boundaryEdge->getLeftFace() == face.first) {
    //             boundaryFaces.insert({dist, boundaryEdge->getRightFace()});
    //         }

    //         else {
    //             boundaryFaces.insert({dist, boundaryEdge->getLeftFace()});
    //         }
    //     }
        
    //     for (auto& targetFace : boundaryFaces) {
    //         if (degFaces[targetFace.second] == 0) {
    //             degFaces[targetFace.second] = 1;
    //             degFaces[face.first] = 1;
                
    //             fout3 << face.first->getGenerator()->getDisk().getX() << "," << face.first->getGenerator()->getDisk().getY() << "," 
    //                   << targetFace.second->getGenerator()->getDisk().getX()  << "," << targetFace.second->getGenerator()->getDisk().getY()  << "\n";
    //             break;
    //         }
    //     }
    // }
    // fout3.close();
    
    // TODO: eulerian path
    // TODO: hamillton path
    
    
}

vector<pair<float, vector<int>>> HeuristicAlgorithm::initialize_chromosome_with_VD(const int& population) {
    vector<pair<float, vector<int>>> chromosomes;
    
    list<VFace2D*> faces; 

    m_VD.getVoronoiFaces(faces);

    faces.pop_front();
    
    for (int i = 0; i < population; ++i) {
        vector<int> tempPath = generate_chromosome(faces);
        chromosomes.push_back({evaluate_function(tempPath), tempPath});
    }
    
    return chromosomes;
}

vector<int> HeuristicAlgorithm::generate_chromosome(const list<VFace2D*>& faces) {
    std::ofstream fout("./../data/result1.txt");
    
    vector<int> result;
    
    bool status = false;
    multimap<int, VFace2D*> chainCountEdges;
    vector<list<VFace2D*>> connectedChain;
    map<VFace2D*, int> connectedFaces;
    
    vector<VFace2D*> faceVector(faces.begin(), faces.end());
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(faceVector.begin(), faceVector.end(), g);
    
    for (const auto& face : faceVector) {
        chainCountEdges.insert({0, face});
    }
    
    auto currFaceIter = chainCountEdges.begin();
    auto currFace = chainCountEdges.begin()->second;
    while (chainCountEdges.count(0) > 0) {        
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
                float dist = currFace->getGenerator()->getDisk().getCenterPt()
                             .distance(candidate->second->getGenerator()->getDisk().getCenterPt());
                tempMap.insert({dist, candidate->second});
            }
            
            auto tempFace = tempMap.begin();
            while (check_same_chain(currFace, tempFace->second, connectedChain) || check_target_face_state(tempFace->second, connectedFaces)) tempFace++;
            
            connect_chain(currFace, tempFace->second, chainCountEdges, connectedChain, connectedFaces);
            fout << currFace->getGenerator()->getDisk().getX() << "," << currFace->getGenerator()->getDisk().getY() << "," 
                << tempFace->second->getGenerator()->getDisk().getX() << "," << tempFace->second->getGenerator()->getDisk().getY() << "\n";

            currFaceIter = chainCountEdges.begin();
            currFace = currFaceIter->second;
        }
        status = false;
        if (chainCountEdges.count(2) > 1 && connectedChain.size() == 1) {
            break;
        }
    }
    for (const auto& city : connectedChain[0]) {
        result.push_back(city->getID()-1);
    }
    fout.close();
    return result;
}

bool HeuristicAlgorithm::check_same_chain(VFace2D* currFace, VFace2D* targetFace, vector<list<VFace2D*>>& connectedChain) {
    for (const auto& currentChain : connectedChain) {
        if ((currentChain.front() == currFace && currentChain.back() == targetFace) ||
            (currentChain.front() == targetFace && currentChain.back() == currFace)) {
                return true;
        }
    }
    return false;
}

bool HeuristicAlgorithm::check_target_face_state(VFace2D* targetFace, map<VFace2D*, int>& connectedFaces) {
    if (connectedFaces[targetFace] == 2) {
        return true;
    }
    return false;
}

void HeuristicAlgorithm::connect_chain(VFace2D* currFace,
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

// improved edge recombination crossover
vector<pair<float, vector<int>>> HeuristicAlgorithm::crossover(vector<pair<float, vector<int>>>& selectionPopulations) { 
    vector<pair<float, vector<int>>*> selectParents = select_parents(selectionPopulations);
    
    while (selectParents.size() > 1) {
        vector<int> gene1 = selectParents[selectParents.size()-2]->second;
        vector<int> gene2 = selectParents[selectParents.size()-1]->second;
        pair<float, vector<int>> offspring1;
        pair<float, vector<int>> offspring2;   

        for (int start = 0; start < 2; ++start) {
            // ------------- make edge
            map<int, vector<int>> edge;
            for (int i = 0; i < gene1.size(); ++i) {
                int previous, next;
                if (i == 0) {
                    previous = gene1.size() - 1;
                    next = i + 1;
                } else if (i == gene1.size() - 1) {
                    previous = i - 1;
                    next = 0;
                } else {
                    previous = i - 1;
                    next = i + 1;
                }

                check_same_value(edge[gene1[i]], gene1[previous]);
                check_same_value(edge[gene1[i]], gene1[next]);
                check_same_value(edge[gene2[i]], gene2[previous]);
                check_same_value(edge[gene2[i]], gene2[next]);
            }
            // ------------- make edge
            vector<int> offspring;
            int currPos = generate_random_int(0, gene1.size()-1);
            erase_value_from_edge(edge, currPos);
            offspring.push_back(currPos);

            for (int i = 0; i < gene1.size() - 1; ++i) {
                multimap<int, int> candidate;
                multimap<int, int> minusCandidate;
           
                for (int n = 0; n < edge[currPos].size(); ++n) {
                    int chromosome = edge[currPos][n];
                    int hasEdgesToNumber = edge[abs(chromosome)].size();

                    if (chromosome >= 0) {
                        candidate.insert({hasEdgesToNumber, chromosome});
                    }

                    else if (chromosome < 0) {
                        minusCandidate.insert({hasEdgesToNumber, chromosome});
                    }
                }

                edge.erase(currPos);

                if (minusCandidate.size() > 0) {
                    if (minusCandidate.size() > 1) {
                        auto it = minusCandidate.begin();
                        if (it->first == (++it)->first) {
                            int num = generate_random_int(0, 1);
                            if (num == 0) {
                                currPos = abs(it->second);
                            } 
                            else {
                                currPos = abs(minusCandidate.begin()->second);
                            }
                        }
                        else {
                            currPos = abs(minusCandidate.begin()->second);
                        }
                    } 
                    else {
                        currPos = abs(minusCandidate.begin()->second);
                    }
                }

                else if (candidate.size() > 0) {
                    if (candidate.size() > 1) {
                        auto it = candidate.begin();
                        if (it->first == (++it)->first) {
                            int num = generate_random_int(0, 1);
                            if (num == 0) {
                                currPos = abs(it->second);
                            } 
                            else {
                                currPos = abs(candidate.begin()->second);
                            }
                        }
                        else {
                             currPos = abs(candidate.begin()->second);
                        }
                    }
                    else {
                        currPos = abs(candidate.begin()->second);
                    }
                }

                else {
                    currPos = edge.begin()->first;
                }
                offspring.push_back(currPos);
                erase_value_from_edge(edge, currPos);
            }
            if (start == 0) {
                offspring1 = {evaluate_function(offspring), offspring};
                mutation(offspring1);
            } 
            else {
                offspring2 = {evaluate_function(offspring), offspring};
                mutation(offspring2);
            }
        }
        *selectParents[selectParents.size() - 2] = offspring1;
        *selectParents[selectParents.size() - 1] = offspring2;
        selectParents.pop_back();
        selectParents.pop_back();
    }
    return selectionPopulations;
}

void HeuristicAlgorithm::order_crossover(vector<pair<float, vector<int>>>& selectionPopulations) {
    vector<pair<float, vector<int>>*> selectParents = select_parents(selectionPopulations);

    while (selectParents.size() > 1) {
        vector<int> gene1 = selectParents[selectParents.size() - 2]->second;
        vector<int> gene2 = selectParents[selectParents.size() - 1]->second;
        pair<float, vector<int>> offspring1;
        pair<float, vector<int>> offspring2;
        
        int start = generate_random_int(0, gene1.size() - 1);
        int end = generate_random_int(start, gene1.size() - 1);
        
        vector<int> temp;
        temp.assign(gene1.begin()+start, gene1.begin()+end);
    }
}


void HeuristicAlgorithm::mutation(pair<float, vector<int>>& offspring) {
    if (m_mutation == "inversion") {
        inversion_mutation(offspring);
    }
    else if (m_mutation == "insertion") {
        insertion_mutation(offspring);
    }
    else if (m_mutation == "displacement") {
        displacement_mutation(offspring);
    }
    else {
        std::cout << "wrong mutation parameter" << std::endl;
        exit(0);
    }
}

void HeuristicAlgorithm::inversion_mutation(pair<float, vector<int>>& crossoverPopulations) {
    float randomProb = generate_random_float(0, 1);
    if (randomProb < m_mutationParameter) {
        int start = generate_random_int(0, crossoverPopulations.second.size() - 1);
        int end = generate_random_int(start, crossoverPopulations.second.size() - 1);
        std::reverse(crossoverPopulations.second.begin() + start, crossoverPopulations.second.begin() + end);
        crossoverPopulations.first = evaluate_function(crossoverPopulations.second);
    }    
}

void HeuristicAlgorithm::insertion_mutation(pair<float, vector<int>>& crossoverPopulations) {
    float randomProb = generate_random_float(0, 1);
    if (randomProb < m_mutationParameter) {
        int start = generate_random_int(0, crossoverPopulations.second.size() - 1);
        int end = start + 1;
        vector<int> temp;
        temp.assign(crossoverPopulations.second.begin()+start, crossoverPopulations.second.begin()+end);
        crossoverPopulations.second.erase(crossoverPopulations.second.begin()+start, crossoverPopulations.second.begin()+end);
        
        int position = generate_random_int(0, crossoverPopulations.second.size() - 1);
        crossoverPopulations.second.insert(crossoverPopulations.second.begin() + position, temp.begin(), temp.end());
        
        crossoverPopulations.first = evaluate_function(crossoverPopulations.second);
    }    
}

void HeuristicAlgorithm::displacement_mutation(pair<float, vector<int>>& crossoverPopulations) {
    float randomProb = generate_random_float(0, 1);
    if (randomProb < m_mutationParameter) {
        int start = generate_random_int(0, crossoverPopulations.second.size() - 1);
        int end = generate_random_int(start, crossoverPopulations.second.size() - 1);
        vector<int> temp;
        temp.assign(crossoverPopulations.second.begin()+start, crossoverPopulations.second.begin()+end);
        crossoverPopulations.second.erase(crossoverPopulations.second.begin()+start, crossoverPopulations.second.begin()+end);
        
        int position = generate_random_int(0, crossoverPopulations.second.size() - 1);
        crossoverPopulations.second.insert(crossoverPopulations.second.begin() + position, temp.begin(), temp.end());
        
        crossoverPopulations.first = evaluate_function(crossoverPopulations.second);
    }   
}

vector<pair<float, vector<int>>> HeuristicAlgorithm::initialize_chromosome(const int& population) {
    vector<pair<float, vector<int>>> chromosomes;
    vector<int> chromosome;
    for (int i = 0; i < m_distanceMatrix.size(); ++i) {
        chromosome.push_back(i);
    }
    
    chromosomes.resize(population);    
    
    for (int i = 0; i < population; ++i) {
        vector<int> tempVector = chromosome;
        for (int j = 0; j < chromosome.size(); ++j) {
            int num = generate_random_int(0, tempVector.size() - 1);
            chromosomes[i].second.push_back(tempVector[num]);
            tempVector.erase(tempVector.begin() + num);
        }
        chromosomes[i].first = evaluate_function(chromosomes[i].second);
    }
    
    return chromosomes;
}

vector<pair<float, vector<int>>> HeuristicAlgorithm::selection(const vector<pair<float, vector<int>>>& chromosomes) {
    vector<pair<float, vector<int>>> selectionChromosome;
    multimap<float, vector<int>> sortedPopulation;
    float maxFitness = 0;
    float bestFitness = 0;
    float worstFitness = 0;
    
    for (auto it = chromosomes.begin(); it != chromosomes.end(); ++it) {
        sortedPopulation.insert({it->first, it->second});
    }
    
    bestFitness = sortedPopulation.begin()->first;
    worstFitness = sortedPopulation.rbegin()->first;
    
    for (auto it = chromosomes.begin(); it != chromosomes.end(); ++it) {
        maxFitness += (worstFitness - it->first) + ((worstFitness - bestFitness) / (m_selectionPressure - 1));
    }
    
    auto it = sortedPopulation.begin();

    int eliteProportion = m_population * m_eliteProportion;
    while(selectionChromosome.size() < eliteProportion && it != sortedPopulation.end()) {
        if (selectionChromosome.size() == 0 || selectionChromosome.rbegin()->first != it->first) {
            selectionChromosome.push_back({it->first, it->second});
        }
        it++;
    }
    int endP = selectionChromosome.size();
    
    for (int i = 0; i < chromosomes.size() - endP; ++i) {
        float randomNumber = generate_random_float(0, maxFitness);
        float sum          = 0;
        
        for (auto it = chromosomes.begin(); it != chromosomes.end(); ++it) {
            float fitness = (worstFitness - it->first) + ((worstFitness - bestFitness) / (m_selectionPressure - 1));

            sum += fitness;

            if (randomNumber < sum) {
                selectionChromosome.push_back({it->first, it->second});
                break;
            }
        }
    }
    
    return selectionChromosome;
}

vector<pair<float, vector<int>>> HeuristicAlgorithm::evaluation(const vector<vector<int>>& populations) {
    vector<pair<float, vector<int>>> result;
    
    for (int i = 0; i < populations.size(); ++i) {
        float fitness = evaluate_function(populations[i]);
        result.push_back({fitness, populations[i]});
    }
    
    return result;
}

float HeuristicAlgorithm::evaluate_function(const vector<int>& population) {
    float result = 0;
    // for (int chromosome = 0; chromosome < population.size(); ++chromosome) {
    //     if (chromosome != population.size() - 1) {
    //         result += m_cities[population[chromosome]].distance_to(m_cities[population[chromosome + 1]]);
    //     }

    //     else {
    //         result += m_cities[population[chromosome]].distance_to(m_cities[population[0]]);
    //     }
    // };
    
    for (int chromosome = 0; chromosome < population.size(); ++chromosome) {
        if (chromosome != population.size() - 1) {
            result += m_distanceMatrix[population[chromosome]][population[chromosome+1]];
        }

        else {
            result += m_distanceMatrix[population[chromosome]][population[0]];
        }
    };

    return result;
}

void HeuristicAlgorithm::check_same_value(vector<int>& edge, const int& value) {
    for (int i = 0; i < edge.size(); ++i) {
        if (abs(edge[i]) == abs(value)) {
            edge[i] = abs(edge[i]) * -1;
            return;
        }
    }
    
    edge.push_back(value);
    return;
}

void HeuristicAlgorithm::erase_value_from_edge(map<int, vector<int>>& edge, const int& value) {
    for (auto it = edge.begin(); it != edge.end(); ++it) {
        for (int i = 0; i < it->second.size(); ++i) {
            if (abs(it->second[i]) == value) {
                it->second.erase(it->second.begin() + i);
            }
        }
    }
}

vector<pair<float, vector<int>>*> HeuristicAlgorithm::select_parents(vector<pair<float, vector<int>>>& selectionPopulations) {
    vector<pair<float, vector<int>>*> result;

    for (auto& chromosome : selectionPopulations) {
        float randomNum = generate_random_float(0, 1);
        if (randomNum <= m_crossoverParameter) {
            result.push_back(&chromosome);
        }
        else if (chromosome.first <= m_bestSolution.first) {
            result.push_back(&chromosome);
        }
    }
    
    return result;
}

pair<float, vector<int>> HeuristicAlgorithm::find_best_fitness(const vector<pair<float, vector<int>>>& populations) {
    for (const auto& it : populations) {
        if (it.first < m_bestSolution.first || m_bestSolution.first == 0) {
            m_bestSolution = it;
            end = clock();
            result = end - start;
            vector<float> info = {m_selectionPressure,
                                  m_crossoverParameter,
                                  m_mutationParameter,
                                  float(m_population),
                                  m_eliteProportion,
                                  float(m_currGeneration),
                                  result / CLOCKS_PER_SEC};
            save_best_solution(info);
        }
    }
    return m_bestSolution;
}

void HeuristicAlgorithm::generate_cities() {
    FileStream file;
    file.read(m_tspFile);
    m_cities = file.get_cities();
    for (int i = 0; i < m_cities.size(); ++i) {
        vector<float> tempVector;
        for (int j = 0; j < m_cities.size(); ++j) {
            tempVector.push_back(m_cities[i].distance_to(m_cities[j]));
        }
        m_distanceMatrix.push_back(tempVector);
    }
}

void HeuristicAlgorithm::generate_distance_matrix() {
    FileStream file;
    file.read_distance_matrix(m_distanceMatrixFile);
    m_distanceMatrix = file.get_distnance_matrix();
}


void HeuristicAlgorithm::save_best_solution(const vector<float>& info) {    
    FileStream file;
    string savePath = m_savePath + m_saveFile;
    file.write(savePath, m_bestSolution, info);
}

int HeuristicAlgorithm::generate_random_int(const int& min, const int& max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(min, max);

    return dis(gen);
}

float HeuristicAlgorithm::generate_random_float(const float& min, const float& max) {
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    
    return dis(gen);
}

void HeuristicAlgorithm::save_vertex(const string& path, const list<VVertex2D*> vertices)  {
    FileStream file;
    file.write_to_vertices(path, vertices);
}
void HeuristicAlgorithm::save_edge(const string& path, const list<VEdge2D*> edges)  {
    FileStream file;
    file.write_to_edge(path, edges);
}
void HeuristicAlgorithm::save_face(const string& path, const list<VFace2D*> faces)  {
    FileStream file;
    file.write_to_face(path, faces);
}