#include "heuristic_algorithm.h"

#include <time.h>

#include <algorithm>
#include <random>
#include <iostream>

#include "file_stream.h"

HeuristicAlgorithm::HeuristicAlgorithm() {
    start = clock();
    generate_cities();
    // generate_distance_matrix();
    generate_vd();
    vector<pair<float, vector<int>>> populations = initialize_chromosome_with_VD(m_population);
    
    // vector<pair<float, vector<int>>> populations = initialize_chromosome(m_population);
    // for (int i = 0; i < m_generation; ++i) {
    //     populations = selection(populations);
    //     populations = crossover(populations);
        
    //     m_currGeneration = i;
    //     if (i % 1000 == 0) {
    //         std::cout << "Generation " << i << " : " << find_best_fitness(populations).first << std::endl;
    //     }
    // }
    // end = clock();
    // result = end - start;
    // vector<float> info = {m_selectionPressure, 
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
    
    // m_VD.constructVoronoiDiagram(circles);
    m_VD.constructVoronoiDiagramCIC_noContainerInInput(circles);
};

vector<pair<float, vector<int>>> HeuristicAlgorithm::initialize_chromosome_with_VD(const int& population) {
    vector<pair<float, vector<int>>> chromosomes;
    
    list<VFace2D*> faces; 
    list<VEdge2D*> edges; 
    list<Generator2D*> generators; 

    m_VD.getGenerators(generators);
    m_VD.getVoronoiFaces(faces);
    m_VD.getVoronoiEdges(edges);

    auto a = m_VD.getContainerGenerator();

    faces.pop_front();
    list<VEdge2D*> resultEdges;
    m_VD.getVoronoiEdges(resultEdges);
    
    for (const auto& face : faces) {
        list<VEdge2D*> edges;
        face->getBoundaryVEdges(edges);
        
        resultEdges.push_back(edges.front());        
    }
    
    FileStream file;
    file.write_to_edge("a", edges);
    file.write_to_face("a", faces);
    // for (int i = 0; i < m_VD.size(); ++i) {
    //     int randomNum = generate_random_int(0, m_cities.size());
        
    //     std::cout << m_
    // }
    
    return chromosomes;
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