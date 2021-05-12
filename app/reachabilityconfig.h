/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Jonathan Trummer
 *      Kathrin Hanauer
 *****************************************************************************/

#pragma once
#include "extern/KaHIP/lib/tools/random_functions.h"
#include <vector>
struct ReachabilityConfig {
    MersenneTwister::result_type seed = time(0);
    unsigned int forward = 5;
    unsigned int binSearchThreshold = 100;

    bool pruning_bibfs = false;
    bool pruning_dfs = false;
    bool edgeCheck = false;
    bool fullReach = false;
    std::string graph_input = "";
    std::string output = "results.csv";
    std::string observations = "observation_dist.csv";
    // std::string inputQueryFile = "";
    std::vector<std::string> queryFiles;
    std::string outputQueryFile = "";
    std::string outputCondensation = "";
    std::string outputUnanswered = "";

    bool useQueryFile = false;
    bool writeOutputQueryFile = false;
    bool writeResultFile = false;
    bool writeObservationsFile = false;
    bool writeCondensationGraph = false;
    bool writeUnansweredQueries = false;

    bool runBfs = false;
    bool runBibfs = false;

    double l = 0.02;
    double l2 = 0.02;
    // support = 1 : sample from topsort range, |R-(v)|*|R+(v)| is maximized
    // support = 2 : n*l nodes sampled, |R-(v)|*|R+(v)| + |R-(v)|*|V\R-(v)| + |R+(v)|*|V\R+(v)| is maximized
    unsigned int support = 3;
    unsigned int t = 8;

};
