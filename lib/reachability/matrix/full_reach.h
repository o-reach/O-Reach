/******************************************************************************
 * This file is part of O'Reach and released under MIT License
 * Authors: 
 *      Jonathan Trummer
 *****************************************************************************/

#pragma once

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/tools/random_functions.h"
#include "app/reachabilityconfig.h"
#include "lib/reachability/Query.h"
#include "lib/algorithms/bfs.h"

#include "extern/flat_hash_map/bytell_hash_map.hpp"
#include "extern/murmur_hash/MurmurHash3.h"
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <omp.h>

struct MurmurHash {
    uint32_t seed{10};
    std::size_t operator()(const NodeID & id) const {
        std::size_t result{0};
        MurmurHash3_x86_32(&id, sizeof(NodeID), seed, &result);
        return result;
    }
};
// using HashFunction = MurmurHash;
using HashFunction = std::hash<NodeID>;
using HashSet = ska::bytell_hash_set<NodeID, HashFunction>;
// using HashSet = std::unordered_set<NodeID>;


class FullReach {
    public: 
        FullReach(graph_access & G): G(G) {
            reach.assign(G.number_of_nodes(), HashSet());
        }

        void query(std::vector<Query> &queries) {
            for (auto & q : queries) {
                q.reachable = reach[q.s].count(q.t);
                q.setResultFromBinary();
            }
        }

        void init() {
            auto num_threads = 1U;
            #pragma omp parallel 
            {
                if (omp_get_thread_num() == 0) {
                    num_threads = omp_get_num_threads();
                }
            }
            std::vector<bfs> bfs_query(num_threads, G);

            #pragma omp parallel for schedule(dynamic)
            for (NodeID u = 0; u < G.number_of_nodes(); u++) {
                auto threadnum = omp_get_thread_num();
                bfs_query[threadnum].explore(u);
                for (auto & n : bfs_query[threadnum].getTouched()) {
                    reach[u].insert(n);
                    // std::cout << n << " " << HashFunction()(n) << "\n";
                }
                bfs_query[threadnum].doCleanup();
            }
            #ifdef RDEBUG
            unsigned int sum{0};
            double load_factors{0};
            unsigned int max_size{0};
            for (auto & n : reach) {
                sum+=n.size();
                load_factors += n.load_factor();
                max_size += n.bucket_count();
            }
            std::cout << "number of NodeIDs stored (ie reach): " << sum << "\n";
            std::cout << "load factor: " << load_factors / reach.size() << "\n";
            std::cout << "max size: " << max_size << "\n";
            #endif
        }

    private:
        graph_access & G;
        
        std::vector<HashSet> reach;

};
