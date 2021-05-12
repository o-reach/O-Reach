/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Jonathan Trummer
 *****************************************************************************/

#pragma once

#include "lib/reachability/Query.h"
#include <vector>
#include "lib/algorithms/bfs.h"
#include "lib/reachability/oreach.h"
#include "app/reachabilityconfig.h"
#include <iostream>
#include "extern/KaHIP/lib/definitions.h"
#include "extern/KaHIP/lib/tools/random_functions.h"
#include <algorithm>


class QueryGenerator {
    ReachabilityConfig & config;
    graph_access & G;
    bfs * bfs_instance;
    unsigned int num_queries;
    bool writeOutput{true};

    std::vector<Query> queries;
    oreach & reach;

    unsigned int n;

    std::vector<char> neighbors;
    std::vector<NodeID> touched;

    inline NodeID randomNode() {
        return random_functions::nextInt(0, n-1);
    }

    inline Query randomQuery(bool preventSameNode = false) {
        NodeID s = randomNode();
        NodeID t = randomNode();
        while (preventSameNode && s==t) {
            t = randomNode();
        }
        return Query(s, t);
    }


    void writeToFile(std::string filename) {
        std::ofstream queryfile;
        queryfile.open(filename, std::ios::out | std::ios::trunc);

        for (auto i = 0U; i < num_queries; i++) {
            queryfile << queries[i].s << " " << queries[i].t << " "
            << (queries[i].reachable ? "1" : "0") << " "
            << (queries[i].answeredWithObservation() ? "1" : "0")
            << "\n";
        }
        queryfile.close();
        std::cout << "wrote queries to " << filename << "\n";
        queries.clear();
    }

    void clearNeighbors() {
        for (auto el : touched) {
            neighbors[el] = false;
        }
        touched.clear();
        for (std::size_t i = 0; i < neighbors.size(); i++) {
            if (neighbors[i]) {
                std::cout << "not reset!! " << i << "\n"; 
            }
        }
    }

    public:
        QueryGenerator(ReachabilityConfig & c, graph_access & G, unsigned int num_queries, oreach & reach): config(c), G(G), num_queries(num_queries), reach(reach) {
            bfs_instance = new bfs(G);
            queries.reserve(num_queries);
            
            n = G.number_of_nodes();
            neighbors.assign(n, false);
        }

        ~QueryGenerator() {
            delete bfs_instance;
        }

        // positive queries using bfs to explore neighbors of a vertex
        void generate_positive() {
            std::cout << "generating " << num_queries << " positive queries \n";
            std::size_t maxEls = (num_queries / n) * 0.75 + 2;
            while (queries.size() < num_queries) {
                NodeID next = randomNode();

                if (G.getNodeDegree(next) == 0) {
                    continue;
                }

                unsigned count = bfs_instance->explore_neighbors(neighbors, next, touched);

                if (count < 1) {
                    clearNeighbors();
                    continue;
                }
                random_functions::permutate_vector_fast(touched, false);
                auto boundary = std::min(maxEls, touched.size());
                for (std::size_t i = 0; i < boundary; i++) {
                    if (touched[i] != next) {
                        Query q(next, touched[i]);
                        q.reachable = true;
                        queries.push_back(q);
                    }
                }
                clearNeighbors();
            }

            if (writeOutput) {
                writeToFile(config.outputQueryFile + ".t.txt");
            }
        }

        // positive queries by creating random pairs of nodes, checking the reachability
        // if reachable => add to queries
        void generate_positive2() {
            std::cout << "generating " << num_queries << " positive queries \n";

            while (queries.size() < num_queries) {
                Query q = randomQuery();
                reach.singleQuery<FallbackMode::FB_PIBFS>(q);
                while (!q.reachable) {
                    q = randomQuery(true);
                    reach.singleQuery<FallbackMode::FB_PIBFS>(q);
                }
                queries.push_back(q);
            }

            if (writeOutput) {
                writeToFile(config.outputQueryFile + ".t.txt");
            }
        }


        void generate_negative() {
            std::cout << "generating " << num_queries << " negative queries \n";
            for (unsigned int i = 0; i < num_queries; ) {
                Query q(randomNode(), randomNode());
                q.reachable = true;
                reach.singleQuery<FallbackMode::FB_PIBFS>(q);
                // prevent endless / extremly long running loops 
                // when q.s reaches a LOT of other vertices
                int count = 25;
                while (q.reachable && (count > 0)) {
                    q.t = randomNode();
                    reach.singleQuery<FallbackMode::FB_PIBFS>(q);
                    count--;
                }
                if (!q.reachable) {
                    queries.push_back(q);
                    i++;
                }
            }

            if (writeOutput) {
                writeToFile(config.outputQueryFile + ".f.txt");
            }
        }

        void generate_random() {
            std::cout << "generating " << num_queries << " random queries \n";
            for (unsigned int i = 0; i < num_queries; i++) {
                Query q(randomNode(), randomNode());
                while (q.s == q.t) {
                    q.t = randomNode();
                }
                queries.push_back(q);
            }

            if (writeOutput) {
                writeToFile(config.outputQueryFile + ".txt");
            }
        }

        void setWriteOutput(bool write) {
            writeOutput = write;
        }

        std::vector<Query> getQueries() {
            return queries;
        }



};
