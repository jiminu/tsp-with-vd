#include <iostream>
#include "heuristic_algorithm.h"

int main(int, char**) {
    string savePath = "./../data/";
    vector<string> largeDataList = {"pr1002",
                                    "pcb3038",
                                    "usa13509",
                                    "d15112",
                                    "rl5934",
                                    "pla33810"};
    vector<string> smallDataList = {"pr144",
                                    "a280",
                                    "pcb442",
                                    "d493",
                                    "vm1084",
                                    "d1291"};
    string algorithm = "2app";
    
    // for (auto str : smallDataList) {
    //     std::cout << str << std::endl;
    //     HeuristicAlgorithm heurisic(savePath, str, algorithm);
    // }
    HeuristicAlgorithm heurisic(savePath, smallDataList[01], algorithm);
    
    std::cout << "Hello, world!\n";
}
