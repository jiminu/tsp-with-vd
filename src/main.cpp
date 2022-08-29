#include <iostream>
#include "heuristic_algorithm.h"

int main(int, char**) {
    string savePath = "./../data/";
    vector<string> largeDataList = {"pr1002",
                                    "pcb3038",
                                    "rl5934",
                                    "usa13509",
                                    "d15112",
                                    "pla33810"};
    vector<string> smallDataList = {"pr144",
                                    "a280",
                                    "pcb442",
                                    "d493",
                                    "vm1084",
                                    "d1291"};
    vector<string> dlist = { "pr144",
                            "a280",
                            "pcb442",
                            "d493",
                            "d1291",
                            "fl1400"};
    
    string algorithm = "2app";
    
    // for (auto str : dlist) {
    //     std::cout << str << std::endl;
    //     HeuristicAlgorithm heurisic(savePath, str, algorithm);
    //     std::cout << std::endl;
    // }
    HeuristicAlgorithm heurisic(savePath, "d1291", algorithm);
    // HeuristicAlgorithm heurisic;
    
    
    std::cout << "Hello, world!\n";
}
